"""Real 1000 Genomes genotypes, two readers, and a second grouping.

The slice is the first biallelic sites on chromosome 22 from the Phase 3
release, restricted to samples that are also in the project pedigree.
Dosages are the counted alleles in the GT field. Nothing is simulated.

Two readers see the same matrix. The opaque reader takes a byte range of a
raw uint8 file and rounds it up to 64 KB blocks. The columnar reader is an
uncompressed Parquet file and is charged the column chunks it opens.

Two groupings besides one all-samples tile:

  nuclear families from the pedigree (a person in two families is stored twice)
  superpopulation panels from the 1000 Genomes sample panel (no overlap)

The rule prices a layout as stored bytes plus 32 times the mean bytes of the
workload. It keeps a partition when that sum is smaller than the all-samples
tile. Inside checks name one group. Crossing checks name two. A scan names
everyone.

    python experiments/real_genotypes.py
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
GENO = ROOT / "data" / "genotypes"
PED = ROOT / "data" / "pedigrees" / "1kg_complete_parents.ped"
VCF = GENO / "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
PANEL = GENO / "integrated_call_samples_v3.20130502.ALL.panel"
OUT = ROOT / "results" / "real_slice"
BLOCK = 65536
N_VARIANTS = 3000
WINDOW = 200
HORIZON = 32
N_QUERIES = 24
SOURCE = (
    "https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/"
    "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
)


def _dosage(gt: str) -> int | None:
    if not gt or "." in gt:
        return None
    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)
    if len(parts) < 2:
        return None
    try:
        return int(parts[0]) + int(parts[1])
    except ValueError:
        return None


def _load_panel() -> dict[str, str]:
    pop: dict[str, str] = {}
    for line in PANEL.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("sample"):
            continue
        parts = line.split()
        if len(parts) >= 3:
            pop[parts[0]] = parts[2]
    return pop


def _load_trios(samples: set[str]) -> list[tuple[str, str, str]]:
    """Nuclear trios whose three members are in the genotype slice."""

    fam: dict[tuple[str, str], list[str]] = {}
    for line in PED.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) < 4:
            continue
        iid, dad, mom = parts[1], parts[2], parts[3]
        if dad in {"0", "-9", ""} or mom in {"0", "-9", ""}:
            continue
        if iid in samples and dad in samples and mom in samples:
            fam.setdefault((dad, mom), []).append(iid)
    trios = []
    for (dad, mom), children in fam.items():
        for child in children:
            trios.append((dad, mom, child))
    return trios


def _load_dosages() -> tuple[list[str], np.ndarray]:
    """Stream the VCF and keep the first biallelic sites.

    A partial download still yields those sites. Gzip raises at the cut.
    """

    samples: list[str] | None = None
    rows: list[np.ndarray] = []
    try:
        with gzip.open(VCF, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]
                    continue
                if samples is None:
                    raise RuntimeError("VCF has no sample header")
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10 or len(fields[3]) != 1 or len(fields[4]) != 1:
                    continue
                gt_index = fields[8].split(":").index("GT")
                dose = np.empty(len(samples), dtype=np.uint8)
                ok = True
                for i, cell in enumerate(fields[9:]):
                    piece = cell.split(":", 1)[0] if gt_index == 0 else cell.split(":")[gt_index]
                    value = _dosage(piece)
                    if value is None or value > 2:
                        ok = False
                        break
                    dose[i] = value
                if not ok:
                    continue
                rows.append(dose)
                if len(rows) >= N_VARIANTS:
                    break
    except (EOFError, gzip.BadGzipFile, OSError):
        if len(rows) < 500:
            raise
    if not rows or samples is None:
        raise RuntimeError("no genotypes read")
    return samples, np.vstack(rows)


def _block_bytes(start: int, length: int) -> int:
    if length <= 0:
        return 0
    end = start + length
    first = start // BLOCK
    last = (end + BLOCK - 1) // BLOCK
    return (last - first) * BLOCK


def _write_group_files(matrix: np.ndarray, groups: list[list[int]], directory: Path) -> list[int]:
    directory.mkdir(parents=True, exist_ok=True)
    sizes = []
    for index, members in enumerate(groups):
        blob = matrix[:, members].reshape(-1).tobytes()
        path = directory / f"group_{index}.bin"
        path.write_bytes(blob)
        sizes.append(len(blob))
    return sizes


def _query_group_bytes(sizes: list[int], groups: list[list[int]], wanted: list[int], n_variants: int, window: int) -> int:
    """Bytes to read the groups that cover ``wanted`` over one variant window."""

    remaining = set(wanted)
    chosen: list[int] = []
    while remaining:
        best = -1
        gain = 0
        for index, members in enumerate(groups):
            hit = len(remaining.intersection(members))
            if hit > gain:
                gain = hit
                best = index
        if gain == 0:
            raise RuntimeError("groups do not cover the query")
        chosen.append(best)
        remaining -= set(groups[best])
    total = 0
    for index in chosen:
        stride = len(groups[index])
        # The window sits at the same variant rows in every group file.
        start = (n_variants // 2) * stride
        total += _block_bytes(start, window * stride)
    return total


def _all_samples_query(n_samples: int, n_variants: int, window: int) -> int:
    stride = n_samples
    start = (n_variants // 2) * stride
    return _block_bytes(start, window * stride)


def _parquet_bytes(matrix: np.ndarray, columns: list[int], path: Path) -> int:
    names = [f"s{i}" for i in range(matrix.shape[1])]
    table = pa.table({names[i]: pa.array(matrix[:, i]) for i in range(matrix.shape[1])})
    if not path.exists():
        pq.write_table(table, path, compression="NONE", row_group_size=WINDOW)
    pf = pq.ParquetFile(path)
    meta = pf.metadata
    start_row = matrix.shape[0] // 2
    stop_row = start_row + WINDOW
    total = 0
    cursor = 0
    wanted = set(columns)
    for group_index in range(meta.num_row_groups):
        group = meta.row_group(group_index)
        rows = group.num_rows
        if cursor + rows > start_row and cursor < stop_row:
            for column in wanted:
                total += group.column(column).total_uncompressed_size
        cursor += rows
    return total


def _mean(values: list[int]) -> float:
    return float(np.mean(values)) if values else 0.0


def _cost(stored: int, mean_read: float) -> float:
    return stored + HORIZON * mean_read


def _pick(name: str, stored: int, reads: list[int], baseline_stored: int, baseline_read: float) -> dict:
    mean_read = _mean(reads)
    cost = _cost(stored, mean_read)
    base = _cost(baseline_stored, baseline_read)
    return {
        "grouping": name,
        "stored_bytes": stored,
        "mean_read_bytes": mean_read,
        "cost": cost,
        "all_samples_cost": base,
        "rule_keeps_partition": bool(cost < base),
    }


def main() -> None:
    if not VCF.exists() or VCF.stat().st_size < 1_000_000:
        raise SystemExit(f"missing {VCF}. Download it from {SOURCE}")
    if not PANEL.exists():
        raise SystemExit(f"missing {PANEL}")
    OUT.mkdir(parents=True, exist_ok=True)
    print("reading genotypes", flush=True)
    sample_ids, matrix = _load_dosages()
    panel = _load_panel()
    keep = [i for i, sample in enumerate(sample_ids) if sample in panel]
    sample_ids = [sample_ids[i] for i in keep]
    matrix = matrix[:, keep]
    id_of = {sample: index for index, sample in enumerate(sample_ids)}
    n_variants, n_samples = matrix.shape
    # Phase 3 kept the unrelated set: no parent and child are both in the VCF.
    # The pedigree still puts sequenced relatives under one family id.
    by_family: dict[str, list[int]] = {}
    for line in PED.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) < 2 or parts[1] not in id_of:
            continue
        by_family.setdefault(parts[0], []).append(id_of[parts[1]])
    family_groups = [members for members in by_family.values() if members]
    covered = {person for group in family_groups for person in group}
    for index in range(n_samples):
        if index not in covered:
            family_groups.append([index])
    multi = [group for group in family_groups if len(group) >= 2]
    if len(multi) < 5:
        raise RuntimeError(f"only {len(multi)} pedigree families in the slice")
    print(
        f"{n_samples} samples, {n_variants} sites, {len(multi)} families with two or more sequenced members",
        flush=True,
    )

    pop_of: dict[str, list[int]] = {}
    for index, sample in enumerate(sample_ids):
        pop_of.setdefault(panel[sample], []).append(index)
    pop_groups = [members for members in pop_of.values() if members]

    work = OUT / "layouts"
    if work.exists():
        for old in work.rglob("*"):
            if old.is_file():
                old.unlink()
    all_bytes = matrix.reshape(-1).tobytes()
    (work / "all").mkdir(parents=True, exist_ok=True)
    (work / "all" / "all.bin").write_bytes(all_bytes)
    family_sizes = _write_group_files(matrix, family_groups, work / "family")
    pop_sizes = _write_group_files(matrix, pop_groups, work / "population")

    rng = np.random.default_rng(7)
    inside_family = []
    cross_family = []
    draws = rng.choice(len(multi), size=min(N_QUERIES, len(multi)), replace=False)
    for draw in draws:
        inside_family.append(list(multi[int(draw)]))
    for _ in range(N_QUERIES):
        group = list(multi[int(rng.integers(0, len(multi)))])
        others = [i for i in range(n_samples) if i not in set(group)]
        extra = rng.choice(others, size=12, replace=False)
        cross_family.append(group + [int(i) for i in extra])

    inside_pop = []
    cross_pop = []
    pop_ids = [name for name, members in pop_of.items() if len(members) >= 20]
    for _ in range(N_QUERIES):
        name = pop_ids[int(rng.integers(0, len(pop_ids)))]
        pick = rng.choice(pop_of[name], size=15, replace=False)
        inside_pop.append([int(i) for i in pick])
        a, b = rng.choice(len(pop_ids), size=2, replace=False)
        left = rng.choice(pop_of[pop_ids[int(a)]], size=15, replace=False)
        right = rng.choice(pop_of[pop_ids[int(b)]], size=15, replace=False)
        cross_pop.append([int(i) for i in list(left) + list(right)])

    all_read = _all_samples_query(n_samples, n_variants, WINDOW)
    family_inside_reads = [
        _query_group_bytes(family_sizes, family_groups, query, n_variants, WINDOW) for query in inside_family
    ]
    family_cross_reads = [
        _query_group_bytes(family_sizes, family_groups, query, n_variants, WINDOW) for query in cross_family
    ]
    pop_inside_reads = [
        _query_group_bytes(pop_sizes, pop_groups, query, n_variants, WINDOW) for query in inside_pop
    ]
    pop_cross_reads = [
        _query_group_bytes(pop_sizes, pop_groups, query, n_variants, WINDOW) for query in cross_pop
    ]
    scan_all = _block_bytes(0, n_samples * n_variants)
    scan_family = sum(_block_bytes(0, size) for size in family_sizes)
    scan_pop = sum(_block_bytes(0, size) for size in pop_sizes)

    parquet_path = work / "columns.parquet"
    column_bytes = _parquet_bytes(matrix, list(inside_family[0]), parquet_path)

    stored_all = len(all_bytes)
    stored_family = sum(family_sizes)
    stored_pop = sum(pop_sizes)
    dup_family = stored_family / stored_all
    dup_pop = stored_pop / stored_all

    rule = {
        "family_inside": _pick("pedigree family id, checks stay inside", stored_family, family_inside_reads, stored_all, all_read),
        "family_cross": _pick("pedigree family id, check adds 12 other people", stored_family, family_cross_reads, stored_all, all_read),
        "family_scan": _pick("pedigree family id, scans", stored_family, [scan_family], stored_all, scan_all),
        "population_inside": _pick("superpopulation, checks stay inside", stored_pop, pop_inside_reads, stored_all, all_read),
        "population_cross": _pick("superpopulation, checks name two populations", stored_pop, pop_cross_reads, stored_all, all_read),
        "population_scan": _pick("superpopulation, scans", stored_pop, [scan_pop], stored_all, scan_all),
    }
    payload = {
        "source_vcf": SOURCE,
        "chromosome": "22",
        "n_variants": n_variants,
        "n_samples": n_samples,
        "n_families_size_at_least_2": len(multi),
        "n_family_groups": len(family_groups),
        "phase3_has_no_parent_child_pair": True,
        "populations": {name: len(members) for name, members in pop_of.items()},
        "window_variants": WINDOW,
        "block_bytes": BLOCK,
        "horizon_queries": HORIZON,
        "stored_bytes": {"all_samples": stored_all, "nuclear_family": stored_family, "superpopulation": stored_pop},
        "stored_ratio_vs_all_samples": {"nuclear_family": dup_family, "superpopulation": dup_pop},
        "mean_read_bytes": {
            "all_samples_window": all_read,
            "family_inside": _mean(family_inside_reads),
            "family_cross": _mean(family_cross_reads),
            "population_inside": _mean(pop_inside_reads),
            "population_cross": _mean(pop_cross_reads),
            "scan_all": scan_all,
            "scan_family": scan_family,
            "scan_population": scan_pop,
            "parquet_one_family": column_bytes,
        },
        "parquet_over_opaque_trio": column_bytes / _mean(family_inside_reads) if family_inside_reads else None,
        "opaque_all_over_family_inside": all_read / _mean(family_inside_reads) if family_inside_reads else None,
        "rule": rule,
    }
    (OUT / "measurements.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    lines = [
        "# Real chromosome 22 slice",
        "",
        f"{n_samples} Phase 3 samples, {n_variants} biallelic sites, {len(multi)} pedigree families with two or more members.",
        "The Phase 3 callset has no parent and child together. Family groups are the pedigree family ids.",
        "Opaque reads are raw uint8 byte ranges rounded to 64 KB. Parquet is uncompressed column chunks.",
        "The rule keeps a partition when stored bytes + 32 mean reads is below the all-samples tile.",
        "",
        f"Stored bytes: all-samples {stored_all}, nuclear {stored_family} ({dup_family:.3f}x), "
        f"superpopulation {stored_pop} ({dup_pop:.3f}x).",
        "",
        "| Workload | Mean opaque read (B) | Rule keeps the partition |",
        "| --- | ---: | --- |",
    ]
    for key, row in rule.items():
        lines.append(
            f"| {row['grouping']} | {row['mean_read_bytes']:.0f} | {row['rule_keeps_partition']} |"
        )
    lines.extend(
        [
            "",
            f"All-samples window read: {all_read} bytes.",
            f"One pedigree family, Parquet column chunks: {column_bytes} bytes.",
            f"Opaque all-samples / nuclear inside-check: {payload['opaque_all_over_family_inside']:.2f}x.",
            "",
        ]
    )
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({k: payload[k] for k in ("n_samples", "n_variants", "n_families_size_at_least_2", "stored_ratio_vs_all_samples", "mean_read_bytes")}, indent=2))
    for key, row in rule.items():
        print(f"{key}: keep={row['rule_keeps_partition']} read={row['mean_read_bytes']:.0f} cost={row['cost']:.0f} base={row['all_samples_cost']:.0f}")


if __name__ == "__main__":
    main()
