"""Hail baseline on the scaled dosage matrices.

Local mode on one machine, not a Hail cluster. The default object is a
locus-ordered MatrixTable, which is Hail's physical design. The family
object is a Hail Table keyed by nuclear-family id, with a person stored
once per family they belong to. Queries use interval filters, so a family
read opens only the groups the same greedy cover would open in the simulator.

    python experiments/hail_baseline.py --config configs/scale.json --out results/hail
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from pedpack.genotype import build_cohort
from pedpack.grouping import nuclear_groups
from pedpack.pedigree import simulate_pedigree
from pedpack.workloads import make_workloads


def _dir_bytes(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _cover(groups: list[list[int]], samples: tuple[int, ...]) -> list[int]:
    remaining = {int(s) for s in samples}
    chosen: list[int] = []
    sets = [set(group) for group in groups]
    while remaining:
        best_i = -1
        best_gain = 0
        for i, members in enumerate(sets):
            gain = len(members & remaining)
            if gain > best_gain:
                best_gain = gain
                best_i = i
        if best_gain == 0:
            break
        chosen.append(best_i)
        remaining -= sets[best_i]
    if remaining:
        raise RuntimeError("nuclear groups do not cover the query")
    return chosen


def _write_tsv(path: Path, genotypes: np.ndarray, ids: list[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(ids))
        handle.write("\n")
        for row in genotypes:
            handle.write("\t".join(str(int(value)) for value in row))
            handle.write("\n")


def _mean(rows: list[dict], key: str) -> float:
    if not rows:
        return 0.0
    return sum(row[key] for row in rows) / len(rows)


def run_seed(hl, seed: int, config: dict, scratch: Path) -> dict:
    rng = np.random.default_rng(seed)
    pedigree = simulate_pedigree(
        n_families=int(config["n_families"]),
        generations=int(config["generations"]),
        mean_children=float(config["mean_children"]),
        cross_rate=float(config["cross_rate"]),
        rng=rng,
        min_people=int(config.get("min_people", 24)),
    )
    cohort = build_cohort(
        pedigree,
        n_variants=int(config["n_variants"]),
        rng=rng,
        bp_spacing=int(config.get("bp_spacing", 100)),
    )
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
    ids = [str(item) for item in cohort.ids.tolist()]
    seed_dir = scratch / f"seed_{seed}"
    if seed_dir.exists():
        shutil.rmtree(seed_dir)
    seed_dir.mkdir(parents=True)

    tsv = seed_dir / "dosages.tsv"
    started = time.perf_counter()
    _write_tsv(tsv, cohort.genotypes, ids)
    tsv_seconds = time.perf_counter() - started
    print(f"  tsv {tsv_seconds:.1f}s", flush=True)

    started = time.perf_counter()
    matrix = hl.import_matrix_table(str(tsv), delimiter="\t", entry_type=hl.tint32, no_header=False)
    matrix = matrix.add_row_index("row_idx")
    matrix = matrix.annotate_rows(vid=hl.int32(matrix.row_idx))
    matrix = matrix.key_rows_by("vid")
    matrix = matrix.key_cols_by(s=matrix.col_id)
    n_row_parts = max(8, cohort.n_variants // 500)
    matrix = matrix.repartition(n_row_parts, shuffle=True)
    default_path = seed_dir / "default.mt"
    matrix.write(str(default_path), overwrite=True)
    default_write = time.perf_counter() - started
    default = hl.read_matrix_table(str(default_path))
    print(f"  default MatrixTable {default_write:.1f}s, {default.n_partitions()} partitions", flush=True)

    family_dir = seed_dir / "families"
    family_dir.mkdir()
    family_paths: list[str] = []
    started = time.perf_counter()
    for index, group in enumerate(groups):
        columns = [ids[sample] for sample in group]
        tsv_path = family_dir / f"family_{index}.tsv"
        _write_tsv(tsv_path, cohort.genotypes[:, list(group)], columns)
        part = hl.import_matrix_table(str(tsv_path), delimiter="\t", entry_type=hl.tint32, no_header=False)
        part = part.add_row_index("row_idx")
        part = part.annotate_rows(vid=hl.int32(part.row_idx))
        part = part.key_rows_by("vid")
        part = part.key_cols_by(s=part.col_id)
        path = str(family_dir / f"family_{index}.mt")
        part.write(path, overwrite=True)
        family_paths.append(path)
        tsv_path.unlink()
        if (index + 1) % 40 == 0:
            print(f"  wrote {index + 1}/{len(groups)} family MatrixTables", flush=True)
    family_write = time.perf_counter() - started
    print(f"  family MatrixTables {family_write:.1f}s, {len(family_paths)} files", flush=True)

    queries = list(workloads["heldout_family"]) + list(workloads["scan"][:1])
    timed = []
    for query in queries:
        if query.workload == "scan":
            chosen = list(range(len(groups)))
            sample_ids = ids
        else:
            chosen = _cover(groups, query.samples)
            sample_ids = [ids[sample] for sample in query.samples]
        wanted = hl.set(set(sample_ids))
        row_interval = [
            hl.Interval(
                start=hl.Struct(vid=int(query.v0)),
                end=hl.Struct(vid=int(query.v1)),
                includes_start=True,
                includes_end=False,
            )
        ]
        started = time.perf_counter()
        subset = hl.filter_intervals(default, row_interval)
        subset = subset.filter_cols(wanted.contains(subset.s))
        hail_entries = int(subset.entries().count())
        hail_seconds = time.perf_counter() - started

        row_only = [
            hl.Interval(
                start=hl.Struct(vid=int(query.v0)),
                end=hl.Struct(vid=int(query.v1)),
                includes_start=True,
                includes_end=False,
            )
        ]
        started = time.perf_counter()
        family_entries = 0
        for group_index in chosen:
            part = hl.read_matrix_table(family_paths[group_index])
            part = hl.filter_intervals(part, row_only)
            family_entries += int(part.entries().count())
        family_seconds = time.perf_counter() - started
        timed.append(
            {
                "name": query.name,
                "workload": query.workload,
                "n_samples": len(query.samples),
                "n_variants": int(query.v1 - query.v0),
                "n_groups": len(chosen),
                "hail_entries": hail_entries,
                "hail_seconds": hail_seconds,
                "family_entries": family_entries,
                "family_seconds": family_seconds,
            }
        )
        print(
            f"  {query.name}: MatrixTable {hail_seconds:.3f}s / {hail_entries} entries, "
            f"family table {family_seconds:.3f}s / {family_entries} entries",
            flush=True,
        )

    duplication = sum(len(group) for group in groups) / cohort.n_samples
    return {
        "seed": seed,
        "n_samples": cohort.n_samples,
        "n_variants": cohort.n_variants,
        "n_nuclear_groups": len(groups),
        "duplication": duplication,
        "tsv_seconds": tsv_seconds,
        "default_write_seconds": default_write,
        "family_write_seconds": family_write,
        "default_bytes": _dir_bytes(default_path),
        "family_bytes": _dir_bytes(family_dir),
        "n_partitions": default.n_partitions(),
        "family_tables": len(family_paths),
        "queries": timed,
        "engine": f"hail {hl.__version__} local[2]",
        "family_object": "one MatrixTable per nuclear family",
    }


def _report(rows: list[dict]) -> str:
    lines = [
        "# Hail baseline, local mode",
        "",
        "Hail 0.2.139 in Spark local mode on one WSL virtual machine. This is",
        "not a multi-node Hail or HDFS cluster. Both objects are MatrixTables",
        "with the same integer entry type. The default table is variant-keyed",
        "and repartitioned along the locus axis. The family side is one",
        "MatrixTable per nuclear family, so a person in two families is stored",
        "twice and a family query opens only the covering tables.",
        "",
        "| Seed | Samples | Variants | Duplication | MT build s | Family build s | MT bytes | Family bytes | Held-out MT s | Held-out family s |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        held = [query for query in row["queries"] if query["workload"] == "heldout_family"]
        lines.append(
            "| {seed} | {n_samples} | {n_variants} | {dup:.3f} | {mt:.1f} | {fam:.1f} | {mtb} | {famb} | {hts:.3f} | {hfs:.3f} |".format(
                seed=row["seed"],
                n_samples=row["n_samples"],
                n_variants=row["n_variants"],
                dup=row["duplication"],
                mt=row["default_write_seconds"],
                fam=row["family_write_seconds"],
                mtb=row["default_bytes"],
                famb=row["family_bytes"],
                hts=_mean(held, "hail_seconds"),
                hfs=_mean(held, "family_seconds"),
            )
        )
    lines.extend(["", "## Queries", ""])
    for row in rows:
        lines.append(f"### Seed {row['seed']}")
        lines.append("")
        lines.append("| Query | Samples | Variants | Groups | MT s | MT entries | Family s | Family entries |")
        lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        for query in row["queries"]:
            lines.append(
                f"| {query['name']} | {query['n_samples']} | {query['n_variants']} | {query['n_groups']} | "
                f"{query['hail_seconds']:.3f} | {query['hail_entries']} | {query['family_seconds']:.3f} | "
                f"{query['family_entries']} |"
            )
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Time a Hail MatrixTable against a family-keyed Hail Table")
    parser.add_argument("--config", type=Path, default=ROOT / "configs" / "scale.json")
    parser.add_argument("--out", type=Path, default=ROOT / "results" / "hail")
    parser.add_argument("--scratch", type=Path, default=Path("/tmp/pedpack-hail"))
    parser.add_argument("--seeds", type=int, nargs="*", default=None)
    args = parser.parse_args()
    try:
        import hail as hl
    except ImportError as exc:
        raise SystemExit("hail is not installed in this interpreter") from exc

    config = json.loads(args.config.read_text(encoding="utf-8"))
    seeds = args.seeds if args.seeds else [int(seed) for seed in config["seeds"]]
    args.scratch.mkdir(parents=True, exist_ok=True)
    hl.init(
        master="local[2]",
        tmp_dir=str(args.scratch / "spark"),
        quiet=True,
        spark_conf={"spark.driver.memory": "3g"},
    )
    print(f"hail {hl.__version__}", flush=True)
    rows = []
    for seed in seeds:
        print(f"seed {seed}", flush=True)
        rows.append(run_seed(hl, seed, config, args.scratch))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "measurements.json").write_text(json.dumps({"seeds": rows}, indent=2), encoding="utf-8")
    (args.out / "report.md").write_text(_report(rows), encoding="utf-8")
    print(f"wrote {args.out / 'report.md'}")


if __name__ == "__main__":
    main()
