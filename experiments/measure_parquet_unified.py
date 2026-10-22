"""One Parquet file for the packed families, against the 12-file pack.

The 12-file layout is the pack the search chose. Column projection already
returns named samples, and that layout was slower because a query opened
8–10 files. This run puts the same duplicated columns into one file, with
the same 800-variant row groups, and into one file whose row group is the
whole variant axis. Seed 7 of configs/scale.json, same workload draw as
experiments/measure_parquet.py.

    python experiments/measure_parquet_unified.py
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


def _family_columns(groups: list[list[int]], chosen: list[int]) -> list[str]:
    names: list[str] = []
    for group_index in chosen:
        for sample in groups[group_index]:
            names.append(column_name(group_index, sample))
    return names


def _time_read(path: Path, columns: list[str], v0: int, v1: int, row_group_size: int) -> tuple[float, int, int, int]:
    started = time.perf_counter()
    compressed, entries, row_groups = read_window(path, columns, v0, v1, row_group_size)
    return time.perf_counter() - started, compressed, entries, row_groups


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
    )
    groups = nuclear_groups(cohort)
    out = ROOT / "results" / "parquet_unified" / f"seed_{seed}"
    out.mkdir(parents=True, exist_ok=True)

    variant_path = out / "variant.parquet"
    started = time.perf_counter()
    write_dosage_file(
        variant_path,
        {str(sample): cohort.genotypes[:, sample] for sample in range(cohort.n_samples)},
        row_group_size,
    )
    variant_write = time.perf_counter() - started

    split_dir = out / "split"
    split_dir.mkdir(exist_ok=True)
    split_files: list[dict] = []
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
        path = split_dir / f"part-{file_index:05d}.parquet"
        write_dosage_file(path, columns, row_group_size)
        split_files.append({"path": path, "groups": group_ids})
    split_write = time.perf_counter() - started

    unified_columns = {
        column_name(group_index, sample): cohort.genotypes[:, sample]
        for group_index, group in enumerate(groups)
        for sample in group
    }
    unified_path = out / "unified.parquet"
    started = time.perf_counter()
    write_dosage_file(unified_path, unified_columns, row_group_size)
    unified_write = time.perf_counter() - started

    wide_path = out / "unified_wide.parquet"
    started = time.perf_counter()
    write_dosage_file(wide_path, unified_columns, cohort.n_variants)
    wide_write = time.perf_counter() - started
    print(
        f"writes variant {variant_write:.2f}s, split {split_write:.2f}s "
        f"({len(split_files)} files), unified {unified_write:.2f}s, wide {wide_write:.2f}s",
        flush=True,
    )

    queries = list(workloads["heldout_family"]) + list(workloads["scan"][:1])
    timed = []
    for query in queries:
        if query.workload == "scan":
            chosen = list(range(len(groups)))
            variant_columns = [str(sample) for sample in range(cohort.n_samples)]
        else:
            chosen = _cover(groups, query.samples)
            variant_columns = [str(sample) for sample in query.samples]
        family_columns = _family_columns(groups, chosen)
        variant_seconds, variant_bytes, variant_entries, variant_groups = _time_read(
            variant_path, variant_columns, query.v0, query.v1, row_group_size
        )

        chosen_set = set(chosen)
        started = time.perf_counter()
        split_bytes = 0
        split_entries = 0
        split_row_groups = 0
        opened = 0
        for packed in split_files:
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
            split_bytes += compressed
            split_entries += entries
            split_row_groups += row_groups
        split_seconds = time.perf_counter() - started

        unified_seconds, unified_bytes, unified_entries, unified_groups = _time_read(
            unified_path, family_columns, query.v0, query.v1, row_group_size
        )
        wide_seconds, wide_bytes, wide_entries, wide_groups = _time_read(
            wide_path, family_columns, query.v0, query.v1, cohort.n_variants
        )
        row = {
            "name": query.name,
            "workload": query.workload,
            "n_groups": len(chosen),
            "variant_seconds": variant_seconds,
            "variant_compressed_bytes": variant_bytes,
            "variant_entries": variant_entries,
            "variant_row_groups": variant_groups,
            "split_seconds": split_seconds,
            "split_compressed_bytes": split_bytes,
            "split_entries": split_entries,
            "split_files": opened,
            "split_row_groups": split_row_groups,
            "unified_seconds": unified_seconds,
            "unified_compressed_bytes": unified_bytes,
            "unified_entries": unified_entries,
            "unified_row_groups": unified_groups,
            "wide_seconds": wide_seconds,
            "wide_compressed_bytes": wide_bytes,
            "wide_entries": wide_entries,
            "wide_row_groups": wide_groups,
        }
        timed.append(row)
        print(
            f"  {query.name}: variant {variant_seconds:.3f}s/{variant_bytes}B, "
            f"split {split_seconds:.3f}s/{split_bytes}B/{opened} files, "
            f"unified {unified_seconds:.3f}s/{unified_bytes}B, "
            f"wide {wide_seconds:.3f}s/{wide_bytes}B",
            flush=True,
        )

    duplication = sum(len(group) for group in groups) / cohort.n_samples
    payload = {
        "seed": seed,
        "n_samples": cohort.n_samples,
        "n_variants": cohort.n_variants,
        "n_groups": len(groups),
        "n_split_files": len(split_files),
        "groups_per_file": groups_per_file,
        "row_group_size": row_group_size,
        "wide_row_group_size": cohort.n_variants,
        "duplication": duplication,
        "compression": "zstd",
        "variant_write_seconds": variant_write,
        "split_write_seconds": split_write,
        "unified_write_seconds": unified_write,
        "wide_write_seconds": wide_write,
        "variant_bytes": variant_path.stat().st_size,
        "split_bytes": sum(item["path"].stat().st_size for item in split_files),
        "unified_bytes": unified_path.stat().st_size,
        "wide_bytes": wide_path.stat().st_size,
        "queries": timed,
    }
    dest = ROOT / "results" / "parquet_unified"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "measurements.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (dest / "report.md").write_text(_report(payload), encoding="utf-8")
    print(f"wrote {dest / 'report.md'}")


def _mean(rows: list[dict], key: str) -> float:
    return sum(row[key] for row in rows) / len(rows)


def _report(row: dict) -> str:
    held = [query for query in row["queries"] if query["workload"] == "heldout_family"]
    lines = [
        "# One Parquet file for the packed families",
        "",
        "Seed 7 of the scaled cohort. Dosages are the same draw as",
        "`experiments/measure_parquet.py`. Compression is Zstandard.",
        "All-samples is one file, one column per person, row groups of",
        f"{row['row_group_size']} variants. Split is the search pack:",
        f"{row['groups_per_file']} nuclear groups per file,",
        f"{row['n_split_files']} files, duplication {row['duplication']:.3f}.",
        "Unified is those same duplicated columns in one file, same row-group",
        "width, so a query opens one footer and still skips row groups outside",
        "its variant window. Wide is that one file with a single row group",
        "covering every variant, so a window cannot skip.",
        "",
        "| Layout | Files | Row group | Build s | Bytes |",
        "| --- | ---: | ---: | ---: | ---: |",
        f"| all-samples | 1 | {row['row_group_size']} | {row['variant_write_seconds']:.2f} | {row['variant_bytes']} |",
        f"| split count:32 | {row['n_split_files']} | {row['row_group_size']} | {row['split_write_seconds']:.2f} | {row['split_bytes']} |",
        f"| unified | 1 | {row['row_group_size']} | {row['unified_write_seconds']:.2f} | {row['unified_bytes']} |",
        f"| unified wide | 1 | {row['wide_row_group_size']} | {row['wide_write_seconds']:.2f} | {row['wide_bytes']} |",
        "",
        "| Query | Groups | Split files | All-samples s | All-samples B | Split s | Split B | Unified s | Unified B | Wide s | Wide B |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for query in row["queries"]:
        lines.append(
            f"| {query['name']} | {query['n_groups']} | {query['split_files']} | "
            f"{query['variant_seconds']:.3f} | {query['variant_compressed_bytes']} | "
            f"{query['split_seconds']:.3f} | {query['split_compressed_bytes']} | "
            f"{query['unified_seconds']:.3f} | {query['unified_compressed_bytes']} | "
            f"{query['wide_seconds']:.3f} | {query['wide_compressed_bytes']} |"
        )
    if held:
        lines.extend(
            [
                "",
                "Mean held-out: "
                f"all-samples {_mean(held, 'variant_seconds'):.3f} s / {_mean(held, 'variant_compressed_bytes'):.0f} B; "
                f"split {_mean(held, 'split_seconds'):.3f} s / {_mean(held, 'split_compressed_bytes'):.0f} B; "
                f"unified {_mean(held, 'unified_seconds'):.3f} s / {_mean(held, 'unified_compressed_bytes'):.0f} B; "
                f"wide {_mean(held, 'wide_seconds'):.3f} s / {_mean(held, 'wide_compressed_bytes'):.0f} B.",
                "",
            ]
        )
    return "\n".join(lines)


if __name__ == "__main__":
    main()
