"""Opaque-block reads on the 3,202-sample chromosome 22 panel.

The VCF is the high-coverage phased panel that put relatives back in.
The pedigree is 1kGP.3202_samples.pedigree_info.txt. Dosages are the counted
alleles in GT. Nothing is simulated.

    python experiments/real_trios.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
VCF = ROOT / "data" / "1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf"
PED = ROOT / "data" / "3202_samples.pedigree_info.txt"
OUT = ROOT / "results" / "real_trios"
N_VARIANTS = 3000
WINDOW = 200
HORIZON = 32
BLOCKS = (65536, 128 * 1024 * 1024)


def trios_in_panel(samples: set[str], pedigree: Path) -> list[tuple[str, str, str]]:
    """Children whose father and mother are both in ``samples``."""

    found: list[tuple[str, str, str]] = []
    for line in pedigree.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("sample"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        child, father, mother = parts[0], parts[1], parts[2]
        if father in {"0", "-9"} or mother in {"0", "-9"}:
            continue
        if child in samples and father in samples and mother in samples:
            found.append((father, mother, child))
    return found


def _dosage(gt: str) -> int | None:
    if not gt or "." in gt:
        return None
    sep = "|" if "|" in gt else "/"
    alleles = gt.split(sep)
    if len(alleles) != 2:
        return None
    try:
        value = int(alleles[0]) + int(alleles[1])
    except ValueError:
        return None
    if value < 0 or value > 2:
        return None
    return value


def load_biallelic(limit: int) -> tuple[list[str], np.ndarray, int]:
    """First ``limit`` biallelic SNVs. Returns samples, dosages, lines scanned."""

    samples: list[str] | None = None
    rows: list[np.ndarray] = []
    scanned = 0
    with VCF.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if samples is None:
                raise RuntimeError("VCF has no sample header")
            scanned += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or len(fields[3]) != 1 or len(fields[4]) != 1 or "," in fields[4]:
                continue
            dose = np.empty(len(samples), dtype=np.uint8)
            ok = True
            for index, cell in enumerate(fields[9:]):
                value = _dosage(cell.split(":", 1)[0])
                if value is None:
                    ok = False
                    break
                dose[index] = value
            if not ok:
                continue
            rows.append(dose)
            if len(rows) >= limit:
                break
    if samples is None or not rows:
        raise RuntimeError("no biallelic sites read")
    return samples, np.vstack(rows), scanned


def _blocks(nbytes: int, block: int) -> int:
    if nbytes <= 0:
        return 0
    return (nbytes + block - 1) // block


def _round(nbytes: int, block: int) -> int:
    return _blocks(nbytes, block) * block


def family_groups(
    sample_ids: list[str], trios: list[tuple[str, str, str]]
) -> list[list[int]]:
    """One group per father-mother pair, plus singletons for everyone else."""

    index_of = {sample: index for index, sample in enumerate(sample_ids)}
    groups: dict[tuple[str, str], list[int]] = {}
    for father, mother, child in trios:
        key = (father, mother)
        group = groups.get(key)
        if group is None:
            group = [index_of[father], index_of[mother]]
            groups[key] = group
        child_index = index_of[child]
        if child_index not in group:
            group.append(child_index)
    covered = {person for group in groups.values() for person in group}
    out = list(groups.values())
    for index in range(len(sample_ids)):
        if index not in covered:
            out.append([index])
    return out


def main() -> None:
    if not VCF.exists():
        raise SystemExit(f"missing {VCF}")
    if not PED.exists():
        raise SystemExit(f"missing {PED}")
    OUT.mkdir(parents=True, exist_ok=True)
    print("reading header", flush=True)
    sample_ids, matrix, scanned = load_biallelic(N_VARIANTS)
    trios = trios_in_panel(set(sample_ids), PED)
    if not trios:
        raise SystemExit("no parent-child pair is fully inside this VCF")
    groups = family_groups(sample_ids, trios)
    multi = [group for group in groups if len(group) >= 3]
    stored_slots = sum(len(group) for group in groups)
    n_variants, n_samples = matrix.shape
    duplication = stored_slots / n_samples
    window = min(WINDOW, n_variants)
    trio = multi[0]
    rng = np.random.default_rng(7)
    outside = [i for i in range(n_samples) if i not in set(trio)]
    outsiders = rng.choice(outside, size=12, replace=False).tolist()
    cross = trio + outsiders

    def cover(wanted: list[int], per_sample: int, block: int) -> int:
        remaining = set(wanted)
        total = 0
        for group in groups:
            hit = remaining.intersection(group)
            if not hit:
                continue
            total += _round(len(group) * per_sample, block)
            remaining -= set(group)
            if not remaining:
                break
        if remaining:
            raise RuntimeError("groups do not cover the query")
        return total

    rows = []
    for block in BLOCKS:
        all_inside = _round(n_samples * window, block)
        family_inside = cover(trio, window, block)
        family_cross = cover(cross, window, block)
        all_scan = _round(n_samples * n_variants, block)
        family_scan = cover(list(range(n_samples)), n_variants, block)
        rows.append(
            {
                "block": block,
                "all_samples_window": all_inside,
                "trio_window": family_inside,
                "trio_plus_12_window": family_cross,
                "all_samples_scan": all_scan,
                "family_files_scan": family_scan,
                "window_ratio_all_over_trio": all_inside / family_inside,
            }
        )
    record = {
        "samples": n_samples,
        "biallelic_sites": n_variants,
        "variant_rows_scanned": scanned,
        "trios": len(trios),
        "nuclear_families": len(multi),
        "duplication": duplication,
        "window": window,
        "first_trio_size": len(trio),
        "blocks": rows,
    }
    (OUT / "measurements.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    lines = [
        "# High-coverage chromosome 22, real trios",
        "",
        f"{n_samples} samples, {n_variants} biallelic SNVs from {scanned} variant rows.",
        f"{len(trios)} parent-child pairs have both parents in the VCF.",
        f"{len(multi)} nuclear families have a sequenced child. Stored slots / samples = {duplication:.3f}.",
        f"The inside check is the first such family ({len(trio)} people) over {window} sites.",
        "",
        "| Block | All-samples window | Trio window | Trio+12 | All-samples scan | Family-file scan | Window ratio |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['block']} | {row['all_samples_window']} | {row['trio_window']} | "
            f"{row['trio_plus_12_window']} | {row['all_samples_scan']} | {row['family_files_scan']} | "
            f"{row['window_ratio_all_over_trio']:.2f} |"
        )
    lines.append("")
    text = "\n".join(lines)
    (OUT / "report.md").write_text(text, encoding="utf-8")
    print(text, flush=True)


if __name__ == "__main__":
    main()
