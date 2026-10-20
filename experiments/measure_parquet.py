"""Parquet measurement of all-samples versus packed nuclear families.

The family side uses the pack the search chose (32 groups per file), not one
Hail MatrixTable per family. Seed 7 of configs/scale.json is the cohort the
361-table Hail run used.

    python experiments/measure_parquet.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from pedpack.columnar import column_name, pack_groups, read_window, write_dosage_file
from pedpack.genotype import build_cohort
from pedpack.grouping import nuclear_groups
from pedpack.pedigree import simulate_pedigree
from pedpack.workloads import make_workloads


def _cover(groups: list[list[int]], samples: tuple[int, ...]) -> list[int]:
    remaining = {int(sample) for sample in samples}
    sets = [set(group) for group in groups]
    chosen: list[int] = []
    while remaining:
        best_i = -1
        best_gain = 0
        for index, members in enumerate(sets):
            gain = len(members & remaining)
            if gain > best_gain:
                best_gain = gain
                best_i = index
        if best_gain == 0:
            break
        chosen.append(best_i)
        remaining -= sets[best_i]
    if remaining:
        raise RuntimeError("nuclear groups do not cover the query")
    return chosen


def _dir_bytes(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def main() -> None:
    config = json.loads((ROOT / "configs" / "scale.json").read_text(encoding="utf-8"))
    seed = 7
    row_group_size = 800
    groups_per_file = 32
    rng = np.random.default_rng(seed)
    pedigree = simulate_pedigree(
        n_families=int(config["n_families"]),
        generations=int(config["generations"]),
        mean_children=float(config["mean_children"]),
        cross_rate=float(config["cross_rate"]),
        rng=rng,
        min_people=int(config.get("min_people", 24)),
    )
    cohort = build_cohort(pedigree, n_variants=int(config["n_variants"]), rng=rng)
    workloads = make_workloads(
        cohort,
        rng,
        n_fit_scan=int(config.get("n_fit_scan", 4)),
        n_fit_family=int(config.get("n_fit_family", 8)),
        n_fit_extended=int(config.get("n_fit_extended", 6)),
        n_eval_scan=int(config.get("n_eval_scan", 4)),
        n_family_selective=int(config.get("n_family_selective", 8)),
        n_extended=int(config.get("n_extended", 6)),
        n_heldout_family=int(config.get("n_heldout_family", 6)),
        n_heldout_extended=int(config.get("n_heldout_extended", 4)),
        family_fraction=config.get("family_fraction"),
    )
    groups = nuclear_groups(cohort)
    out = ROOT / "results" / "parquet" / f"seed_{seed}"
    if out.exists():
        import shutil

        shutil.rmtree(out)
    out.mkdir(parents=True)
    variant_path = out / "variant.parquet"
    started = time.perf_counter()
    write_dosage_file(
        variant_path,
        {str(sample): cohort.genotypes[:, sample] for sample in range(cohort.n_samples)},
        row_group_size,
    )
    variant_write = time.perf_counter() - started

    files: list[dict] = []
    family_dir = out / "family"
    family_dir.mkdir()
    started = time.perf_counter()
    for file_index, packed in enumerate(pack_groups(groups, groups_per_file)):
        columns: dict[str, np.ndarray] = {}
        group_ids: list[int] = []
        base = file_index * groups_per_file
        for offset, group in enumerate(packed):
            group_index = base + offset
            group_ids.append(group_index)
            for sample in group:
                columns[column_name(group_index, sample)] = cohort.genotypes[:, sample]
        path = family_dir / f"part-{file_index:05d}.parquet"
        write_dosage_file(path, columns, row_group_size)
        files.append({"path": path, "groups": group_ids})
    family_write = time.perf_counter() - started
    print(f"variant write {variant_write:.2f}s, family write {family_write:.2f}s, files {len(files)}", flush=True)

    queries = list(workloads["heldout_family"]) + list(workloads["scan"][:1])
    timed = []
    for query in queries:
        if query.workload == "scan":
            chosen = list(range(len(groups)))
            variant_columns = [str(sample) for sample in range(cohort.n_samples)]
        else:
            chosen = _cover(groups, query.samples)
            variant_columns = [str(sample) for sample in query.samples]
        started = time.perf_counter()
        variant_bytes, variant_entries, variant_groups = read_window(
            variant_path, variant_columns, query.v0, query.v1, row_group_size
        )
        variant_seconds = time.perf_counter() - started

        chosen_set = set(chosen)
        started = time.perf_counter()
        family_bytes = 0
        family_entries = 0
        family_row_groups = 0
        opened = 0
        for packed in files:
            names = [
                column_name(group_index, sample)
                for group_index in packed["groups"]
                if group_index in chosen_set
                for sample in groups[group_index]
            ]
            if not names:
                continue
            opened += 1
            compressed, entries, row_groups = read_window(
                packed["path"], names, query.v0, query.v1, row_group_size
            )
            family_bytes += compressed
            family_entries += entries
            family_row_groups += row_groups
        family_seconds = time.perf_counter() - started
        timed.append(
            {
                "name": query.name,
                "workload": query.workload,
                "n_samples": len(query.samples),
                "n_variants": int(query.v1 - query.v0),
                "n_groups": len(chosen),
                "variant_seconds": variant_seconds,
                "variant_compressed_bytes": variant_bytes,
                "variant_entries": variant_entries,
                "variant_row_groups": variant_groups,
                "family_seconds": family_seconds,
                "family_compressed_bytes": family_bytes,
                "family_entries": family_entries,
                "family_files": opened,
                "family_row_groups": family_row_groups,
            }
        )
        print(
            f"  {query.name}: variant {variant_seconds:.3f}s / {variant_bytes} B, "
            f"packed {family_seconds:.3f}s / {family_bytes} B / {opened} files",
            flush=True,
        )

    duplication = sum(len(group) for group in groups) / cohort.n_samples
    payload = {
        "seed": seed,
        "n_samples": cohort.n_samples,
        "n_variants": cohort.n_variants,
        "n_groups": len(groups),
        "n_files": len(files),
        "groups_per_file": groups_per_file,
        "row_group_size": row_group_size,
        "duplication": duplication,
        "compression": "zstd",
        "variant_write_seconds": variant_write,
        "family_write_seconds": family_write,
        "variant_bytes": variant_path.stat().st_size,
        "family_bytes": _dir_bytes(family_dir),
        "queries": timed,
        "note": "Packed files match count:32. This is not one MatrixTable per family.",
    }
    dest = ROOT / "results" / "parquet"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "measurements.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (dest / "report.md").write_text(_report(payload), encoding="utf-8")
    print(f"wrote {dest / 'report.md'}")


def _report(row: dict) -> str:
    held = [query for query in row["queries"] if query["workload"] == "heldout_family"]
    lines = [
        "# Packed Parquet, one host",
        "",
        "Seed 7 of the scaled cohort. Both layouts are Zstandard Parquet with",
        f"row groups of {row['row_group_size']} variants. The all-samples file",
        "stores every person once. The family side packs nuclear groups",
        f"{row['groups_per_file']} to a file ({row['n_files']} files, duplication",
        f"{row['duplication']:.3f}). A query reads the covering groups' columns",
        "and the row groups that overlap its variant window. This is the pack",
        "the search chose. It is not one Hail MatrixTable per family.",
        "",
        "| Layout | Build s | Bytes |",
        "| --- | ---: | ---: |",
        f"| all-samples | {row['variant_write_seconds']:.2f} | {row['variant_bytes']} |",
        f"| packed families | {row['family_write_seconds']:.2f} | {row['family_bytes']} |",
        "",
        "| Query | Groups | Files | Variant s | Variant compressed B | Packed s | Packed compressed B |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for query in row["queries"]:
        lines.append(
            f"| {query['name']} | {query['n_groups']} | {query['family_files']} | "
            f"{query['variant_seconds']:.3f} | {query['variant_compressed_bytes']} | "
            f"{query['family_seconds']:.3f} | {query['family_compressed_bytes']} |"
        )
    if held:
        lines.extend(
            [
                "",
                f"Mean held-out seconds: all-samples {sum(q['variant_seconds'] for q in held) / len(held):.3f}, "
                f"packed {sum(q['family_seconds'] for q in held) / len(held):.3f}.",
                f"Mean held-out compressed bytes: all-samples "
                f"{sum(q['variant_compressed_bytes'] for q in held) / len(held):.0f}, "
                f"packed {sum(q['family_compressed_bytes'] for q in held) / len(held):.0f}.",
                "",
            ]
        )
    return "\n".join(lines)


if __name__ == "__main__":
    main()
