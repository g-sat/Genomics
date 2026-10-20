"""Duplication against family-query fraction.

Writes the closed-form surface in ``pedpack.boundary`` and a few simulator
checks of the same two layouts (all-samples versus nuclear families).

    python experiments/boundary_grid.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib.pyplot as plt
import numpy as np

from pedpack.boundary import family_versus_variant
from pedpack.cost import storage_weight
from pedpack.genotype import build_cohort
from pedpack.grouping import nuclear_groups
from pedpack.layout import prepare_layout
from pedpack.pedigree import simulate_pedigree
from pedpack.selector import measured_objective
from pedpack.simulator import measure_many
from pedpack.study import _fit_weights, make_cluster
from pedpack.types import Candidate
from pedpack.workloads import make_workloads

FRACTIONS = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60]
DUPLICATIONS = [1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50]
CHECKS = (
    {"name": "two_gen", "n_families": 24, "generations": 2, "cross_rate": 0.25, "min_people": 50},
    {"name": "three_gen", "n_families": 16, "generations": 3, "cross_rate": 0.6, "min_people": 50},
)


def _candidate(family: str, groups: list[list[int]], width: int) -> Candidate:
    return Candidate(
        label=family,
        family=family,
        groups=tuple(tuple(group) for group in groups),
        width=width,
        pack="count:32",
        order="group_major",
    )


def _measure_pair(seed: int, spec: dict, fraction: float, cluster, config: dict) -> dict:
    rng = np.random.default_rng(seed)
    pedigree = simulate_pedigree(
        n_families=spec["n_families"],
        generations=spec["generations"],
        mean_children=2.2,
        cross_rate=spec["cross_rate"],
        rng=rng,
        min_people=spec["min_people"],
    )
    cohort = build_cohort(pedigree, n_variants=1600, rng=rng)
    groups = nuclear_groups(cohort)
    duplication = sum(len(group) for group in groups) / cohort.n_samples
    variant = prepare_layout(_candidate("variant", [list(range(cohort.n_samples))], 1600), cohort, cluster)
    family = prepare_layout(_candidate("family", groups, 1600), cohort, cluster)
    local = dict(config)
    local["family_fraction"] = fraction
    workloads = make_workloads(
        cohort,
        rng,
        n_fit_scan=2,
        n_fit_family=4,
        n_fit_extended=2,
        n_eval_scan=1,
        n_family_selective=2,
        n_extended=1,
        n_heldout_family=2,
        n_heldout_extended=1,
        family_fraction=fraction,
    )
    weights = _fit_weights(workloads["fit"], local)
    scan = workloads["scan"][0]
    reference = measure_many(variant, [scan], cluster)[0].seconds
    weight = storage_weight(reference, variant.stats.storage_bytes, cluster.storage_scan_equivalents)
    variant_cost = measured_objective(variant, measure_many(variant, workloads["fit"], cluster), weight, weights)
    family_cost = measured_objective(family, measure_many(family, workloads["fit"], cluster), weight, weights)
    return {
        "seed": seed,
        "pedigree": spec["name"],
        "n_samples": cohort.n_samples,
        "duplication": duplication,
        "family_fraction": fraction,
        "variant_objective": variant_cost,
        "family_objective": family_cost,
        "ratio_variant_over_family": variant_cost / family_cost,
        "family_wins": family_cost < variant_cost,
    }


def main() -> None:
    config = json.loads((ROOT / "configs" / "scale.json").read_text(encoding="utf-8"))
    cluster = make_cluster(config["cluster"])
    surface = []
    for duplication in DUPLICATIONS:
        for fraction in FRACTIONS:
            row = family_versus_variant(1200, 8000, duplication, fraction, cluster)
            surface.append(
                {
                    "duplication": duplication,
                    "family_fraction": fraction,
                    "ratio": row["ratio_variant_over_family"],
                    "family_wins": bool(row["family_wins"]),
                }
            )
    checks = []
    for spec in CHECKS:
        for fraction in (0.05, 0.15, 0.40):
            for seed in (7, 11):
                print(f"measure {spec['name']} f={fraction} seed={seed}", flush=True)
                checks.append(_measure_pair(seed, spec, fraction, cluster, config))
    out = ROOT / "results" / "boundary"
    out.mkdir(parents=True, exist_ok=True)
    (out / "grid.json").write_text(json.dumps({"surface": surface, "checks": checks}, indent=2), encoding="utf-8")
    _plot(surface, checks, out / "heatmap.png")
    _report(surface, checks, out / "report.md")
    print(f"wrote {out / 'report.md'}")


def _plot(surface: list[dict], checks: list[dict], path: Path) -> None:
    fractions = FRACTIONS
    dups = DUPLICATIONS
    grid = np.zeros((len(dups), len(fractions)))
    for row in surface:
        i = dups.index(row["duplication"])
        j = fractions.index(row["family_fraction"])
        grid[i, j] = row["ratio"]
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    image = ax.imshow(grid, origin="lower", aspect="auto", cmap="coolwarm_r", vmin=0.85, vmax=1.15)
    ax.set_xticks(range(len(fractions)), [f"{100 * f:.0f}%" for f in fractions])
    ax.set_yticks(range(len(dups)), [f"{d:.2f}" for d in dups])
    ax.set_xlabel("Family-query fraction")
    ax.set_ylabel("Nuclear duplication")
    ax.set_title("Variant objective / family objective")
    fig.colorbar(image, ax=ax, label="Ratio (above 1: family wins)")
    for check in checks:
        # nearest cell, marked only as a point in data coordinates of the axes
        j = min(range(len(fractions)), key=lambda k: abs(fractions[k] - check["family_fraction"]))
        i = min(range(len(dups)), key=lambda k: abs(dups[k] - check["duplication"]))
        ax.scatter(j, i, s=28, c="black" if check["family_wins"] else "white", edgecolors="black", linewidths=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def _report(surface: list[dict], checks: list[dict], path: Path) -> None:
    wins = [row for row in surface if row["family_wins"]]
    lines = [
        "# Duplication against selectivity",
        "",
        "The surface is `family_versus_variant` at 1,200 samples and 8,000 variants,",
        "using the scaled cluster (400 bytes per genotype, 64 KB blocks, replication 3).",
        "A ratio above 1 means the all-samples layout costs more. Black points on",
        "`heatmap.png` are simulator checks where the nuclear layout won; white",
        "points are checks where all-samples won.",
        "",
        f"Family wins on {len(wins)} of {len(surface)} formula cells.",
        "",
        "## Simulator checks",
        "",
        "| Pedigree | Seed | People | Duplication | Fraction | Variant / family | Family wins |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in checks:
        lines.append(
            f"| {row['pedigree']} | {row['seed']} | {row['n_samples']} | {row['duplication']:.3f} | "
            f"{row['family_fraction']:.2f} | {row['ratio_variant_over_family']:.3f} | {row['family_wins']} |"
        )
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
