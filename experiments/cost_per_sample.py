"""Partitioning versus cost per sample for a hereditary check.

A hereditary check names a child and both parents inside a variant window.
Two partitions of the same dosage matrix are priced per person in the cohort:

  storage seconds per sample = storage term of the fit objective / N
  query seconds per sample   = weighted fit-query seconds / N

The family partition is the nuclear-family oracle. The all-samples partition
is the variant oracle. Both numbers are taken from results already on disk.
The right-hand panel is the shuffle the fit mix pays when a check crosses
the partition. That shuffle is why a byte-only payback, which counts only
the held-out family read, would keep family tiles on the three-generation
cohort where the timed cost per sample goes up.

  python experiments/cost_per_sample.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "cost_per_sample"

STUDIES = (
    ("58–76", ROOT / "results" / "pedigree" / "study.json"),
    ("~1,200", ROOT / "results" / "scale" / "study.json"),
    ("~3,300", ROOT / "results" / "scale_3k" / "study.json"),
    ("3-gen", ROOT / "results" / "scale_overlap" / "study.json"),
    ("1kg ped.", ROOT / "results" / "public_1kg" / "study.json"),
)


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _seconds(seed: dict, fit: dict, horizon: float) -> tuple[float, float]:
    """Storage and per-query seconds. Older reports stored only the bytes."""

    if "storage_seconds" in fit and "weighted_query_seconds" in fit:
        return float(fit["storage_seconds"]), float(fit["weighted_query_seconds"])
    storage = float(fit["storage_bytes"]) * float(seed["storage_seconds_per_byte"])
    query = (float(fit["objective"]) - float(fit["build_seconds"]) - storage) / horizon
    return storage, query


def _rows() -> list[dict]:
    rows: list[dict] = []
    for cohort, path in STUDIES:
        report = _load(path)
        horizon = float(report["config"].get("horizon_queries", 32))
        for seed in report["seeds"]:
            n = int(seed["n_samples"])
            variant = seed["methods"]["variant_oracle"]["workloads"]["fit"]
            family = seed["methods"]["family_oracle"]["workloads"]["fit"]
            held_v = seed["methods"]["variant_oracle"]["workloads"]["heldout_family"]
            held_f = seed["methods"]["family_oracle"]["workloads"]["heldout_family"]
            choice = seed["methods"]["pedpack"]["workloads"]["fit"]["family"]
            variant_storage, variant_query = _seconds(seed, variant, horizon)
            family_storage, family_query = _seconds(seed, family, horizon)
            rows.append(
                {
                    "cohort": cohort,
                    "seed": int(seed["seed"]),
                    "n_samples": n,
                    "n_trios": int(seed["n_trios"]),
                    "pedpack_choice": choice,
                    "family_duplication": float(family["duplication"]),
                    "storage_ratio": (family_storage / n) / (variant_storage / n),
                    "query_ratio": (family_query / n) / (variant_query / n),
                    "shuffle_bytes_per_sample": float(family["mean_shuffle"]) / n,
                    "heldout_byte_ratio": float(held_f["mean_bytes"]) / float(held_v["mean_bytes"]),
                    "stored_byte_ratio": float(family["primary_bytes"]) / float(variant["primary_bytes"]),
                }
            )
    return rows


def _hdfs(n_by_seed: dict[int, int]) -> list[dict]:
    measurements = _load(ROOT / "results" / "hdfs" / "measurements.json")
    rows = []
    for row in measurements:
        if row["scenario"] != "select_05":
            continue
        if row["method"] not in ("pedpack", "pedpack_no_pedigree_features"):
            continue
        seed = int(str(row["seed_dir"]).split("_")[1])
        n = n_by_seed[seed]
        rows.append(
            {
                "seed": seed,
                "n_samples": n,
                "partition": "nuclear family" if row["family"] == "family" else "all-samples",
                "primary_per_sample": float(row["primary_bytes"]) / n,
                "space_per_sample": float(row["space_consumed_bytes"]) / n,
                "check_bytes_per_sample": float(row["heldout_family_mean_bytes"]) / n,
            }
        )
    return rows


def _write(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _median(rows: list[dict], cohort: str, key: str) -> float:
    return float(np.median([row[key] for row in rows if row["cohort"] == cohort]))


def _plot(rows: list[dict], path: Path) -> None:
    cohorts = []
    for row in rows:
        if row["cohort"] not in cohorts:
            cohorts.append(row["cohort"])
    x = np.arange(len(cohorts))
    storage = [_median(rows, cohort, "storage_ratio") for cohort in cohorts]
    query = [_median(rows, cohort, "query_ratio") for cohort in cohorts]
    shuffle = [_median(rows, cohort, "shuffle_bytes_per_sample") / 1e3 for cohort in cohorts]
    raw_choice = [next(row["pedpack_choice"] for row in rows if row["cohort"] == cohort) for cohort in cohorts]
    choice = ["all-samples" if name == "variant" else "family" for name in raw_choice]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.3))
    width = 0.36
    axes[0].bar(x - width / 2, storage, width, color="#4c78a8", label="Storage seconds / sample")
    axes[0].bar(x + width / 2, query, width, color="#f58518", label="Fit-query seconds / sample")
    axes[0].axhline(1.0, color="black", linewidth=0.8)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(cohorts)
    axes[0].set_ylabel("Nuclear family / all-samples")
    axes[0].set_ylim(0, 1.7)
    axes[0].set_title("Cost per sample of each partition")
    axes[0].legend(frameon=False, fontsize=8)
    for i, name in enumerate(choice):
        axes[0].text(i, 1.55, name, ha="center", va="bottom", fontsize=8, color="#333333")
    colors = ["#4c78a8" if name == "family" else "#bab0ac" for name in raw_choice]
    axes[1].bar(x, shuffle, color=colors)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(cohorts)
    axes[1].set_ylabel("Shuffle on the fit mix (KB / sample)")
    axes[1].set_title("Hereditary check that crosses the partition")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = _rows()
    pedigree = _load(ROOT / "results" / "pedigree" / "study.json")
    n_by_seed = {int(seed["seed"]): int(seed["n_samples"]) for seed in pedigree["seeds"]}
    hdfs = _hdfs(n_by_seed)
    _write(OUT / "ratios.csv", rows)
    _write(OUT / "hdfs_per_sample.csv", hdfs)
    figure = OUT / "cost_per_sample.png"
    _plot(rows, figure)
    (ROOT / "paper" / "figures" / "cost_per_sample.png").write_bytes(figure.read_bytes())
    for row in rows:
        print(
            f"{row['cohort']} seed {row['seed']}: storage {row['storage_ratio']:.3f}x, "
            f"query {row['query_ratio']:.3f}x, shuffle {row['shuffle_bytes_per_sample']:.0f} B/sample, "
            f"held-out bytes {row['heldout_byte_ratio']:.3f}x, PedPack {row['pedpack_choice']}"
        )
    print("--- hdfs select_05, bytes per sample ---")
    for row in hdfs:
        print(
            f"seed {row['seed']} {row['partition']}: stored {row['primary_per_sample']:.0f}, "
            f"replicated {row['space_per_sample']:.0f}, check {row['check_bytes_per_sample']:.0f}"
        )


if __name__ == "__main__":
    main()
