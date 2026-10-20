"""Text and figure summaries of a study run."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

WORKLOADS = ("fit", "scan", "selective", "heldout", "window")
MAIN_METHODS = (
    "heuristic_fixed",
    "variant_oracle",
    "family_oracle",
    "generic_oracle",
    "pedigree_oracle",
    "full_oracle",
    "random_budget",
    "analytical_budget",
    "pedpack_generic_pool",
    "pedpack_no_pedigree_features",
    "pedpack",
    "pedpack_tuned_on_scan",
    "pedpack_tuned_on_selective",
)


def _median(values: list[float]) -> float:
    return float(np.median(np.asarray(values, dtype=np.float64)))


def aggregate(report: dict) -> dict:
    seeds = report["seeds"]
    methods = [name for name in MAIN_METHODS if name in seeds[0]["methods"]]
    summary: dict = {"workloads": {}, "model": {}}
    for workload in WORKLOADS:
        summary["workloads"][workload] = {}
        for method in methods:
            costs = [seed["methods"][method]["workloads"][workload]["objective"] for seed in seeds]
            generic = [
                seed["methods"]["generic_oracle"]["workloads"][workload]["objective"] for seed in seeds
            ]
            full = [seed["methods"]["full_oracle"]["workloads"][workload]["objective"] for seed in seeds]
            seconds = [seed["methods"][method]["workloads"][workload]["mean_seconds"] for seed in seeds]
            read = [seed["methods"][method]["workloads"][workload]["mean_bytes"] for seed in seeds]
            remote = [seed["methods"][method]["workloads"][workload]["remote_fraction"] for seed in seeds]
            shuffle = [seed["methods"][method]["workloads"][workload]["mean_shuffle"] for seed in seeds]
            summary["workloads"][workload][method] = {
                "objective_median": _median(costs),
                "objective_min": float(min(costs)),
                "objective_max": float(max(costs)),
                "ratio_to_generic_median": _median([c / g for c, g in zip(costs, generic)]),
                "ratio_to_full_oracle_median": _median([c / f for c, f in zip(costs, full)]),
                "mean_seconds_median": _median(seconds),
                "mean_bytes_median": _median(read),
                "remote_fraction_median": _median(remote),
                "mean_shuffle_median": _median(shuffle),
            }
    mape_n = [seed["model"]["seconds_mape_neural"] for seed in seeds]
    mape_a = [seed["model"]["seconds_mape_analytical"] for seed in seeds]
    spear_n = [seed["model"]["spearman_neural"] for seed in seeds]
    spear_a = [seed["model"]["spearman_analytical"] for seed in seeds]
    summary["model"] = {
        "seconds_mape_neural_median": _median(mape_n),
        "seconds_mape_analytical_median": _median(mape_a),
        "spearman_neural_median": _median(spear_n),
        "spearman_analytical_median": _median(spear_a),
    }
    return summary


def _fmt(value: float, digits: int = 3) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "nan"
    if abs(value) >= 1000 or (abs(value) > 0 and abs(value) < 0.01):
        return f"{value:.3e}"
    return f"{value:.{digits}f}"


def to_markdown(report: dict) -> str:
    summary = aggregate(report)
    lines = [
        "# Simulator report",
        "",
        "Numbers are from the block simulator in this repository. They are not",
        "HDFS or Hail measurements. A ratio below 1 means a lower objective than",
        "the reference on that row. The fit column amortizes build and storage",
        "over a stated mix of scans and family-scoped queries. Other columns add",
        "the full build to that workload's own query times.",
        "",
        f"Seeds: {', '.join(str(seed['seed']) for seed in report['seeds'])}.",
        "",
        "## Objective ratios",
        "",
        "Median across seeds. `generic` is the best index-chunk layout after an",
        "exhaustive search on the fit mix. `oracle` is the best layout in the",
        "whole pool on the fit mix, then frozen.",
        "",
        "| Method | Fit / generic | Scan / generic | Selective / generic | Held-out / generic | Window / generic | Fit / oracle |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for method in summary["workloads"]["fit"]:
        fit = summary["workloads"]["fit"][method]
        lines.append(
            "| {method} | {fit} | {scan} | {sel} | {hold} | {win} | {oracle} |".format(
                method=method,
                fit=_fmt(fit["ratio_to_generic_median"]),
                scan=_fmt(summary["workloads"]["scan"][method]["ratio_to_generic_median"]),
                sel=_fmt(summary["workloads"]["selective"][method]["ratio_to_generic_median"]),
                hold=_fmt(summary["workloads"]["heldout"][method]["ratio_to_generic_median"]),
                win=_fmt(summary["workloads"]["window"][method]["ratio_to_generic_median"]),
                oracle=_fmt(fit["ratio_to_full_oracle_median"]),
            )
        )
    lines.extend(
        [
            "",
            "## Fit-mix execution (median)",
            "",
            "| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for method, row in summary["workloads"]["fit"].items():
        lines.append(
            f"| {method} | {_fmt(row['objective_median'])} | {_fmt(row['mean_seconds_median'])} | "
            f"{_fmt(row['mean_bytes_median'])} | {_fmt(row['remote_fraction_median'])} | "
            f"{_fmt(row['mean_shuffle_median'])} |"
        )
    model = summary["model"]
    lines.extend(
        [
            "",
            "## Cost model, layouts held out of the pilot",
            "",
            f"- Neural seconds MAPE: {_fmt(model['seconds_mape_neural_median'])}",
            f"- Analytical seconds MAPE: {_fmt(model['seconds_mape_analytical_median'])}",
            f"- Neural objective Spearman: {_fmt(model['spearman_neural_median'])}",
            f"- Analytical objective Spearman: {_fmt(model['spearman_analytical_median'])}",
            "",
            "## Layouts chosen on the fit mix",
            "",
        ]
    )
    for seed in report["seeds"]:
        choice = seed["methods"]["pedpack"]["choice"]
        generic = seed["methods"]["generic_oracle"]["choice"]
        lines.append(
            f"- Seed {seed['seed']}: {seed['n_samples']} samples, {seed['n_trios']} trios, "
            f"{seed['n_layouts']} layouts. PedPack chose `{choice['label']}` "
            f"(trio co-location {choice['trio_colocation']:.2f}, {choice['files']} files). "
            f"Generic oracle chose `{generic['label']}` "
            f"(trio co-location {generic['trio_colocation']:.2f}, {generic['files']} files)."
        )
    lines.extend(["", "## Per-seed fit objectives", ""])
    header = "| Seed | " + " | ".join(summary["workloads"]["fit"]) + " |"
    sep = "| --- | " + " | ".join("---:" for _ in summary["workloads"]["fit"]) + " |"
    lines.extend([header, sep])
    for seed in report["seeds"]:
        cells = [
            _fmt(seed["methods"][method]["workloads"]["fit"]["objective"])
            for method in summary["workloads"]["fit"]
        ]
        lines.append(f"| {seed['seed']} | " + " | ".join(cells) + " |")
    lines.extend(
        [
            "",
            "Held-out and window queries were not used to choose `pedpack`.",
            "`best_on_heldout` in the JSON is an analysis upper bound that did see them.",
            "",
        ]
    )
    return "\n".join(lines) + "\n"


def write_report(report: dict, directory: str | Path) -> None:
    path = Path(directory)
    path.mkdir(parents=True, exist_ok=True)
    payload = dict(report)
    payload["aggregate"] = aggregate(report)
    (path / "study.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (path / "report.md").write_text(to_markdown(payload), encoding="utf-8")
    _write_figure(payload, path / "fit_ratios.png")


def _write_figure(report: dict, path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    summary = report["aggregate"]
    methods = [
        name
        for name in (
            "variant_oracle",
            "family_oracle",
            "generic_oracle",
            "heuristic_fixed",
            "analytical_budget",
            "pedpack_generic_pool",
            "pedpack",
        )
        if name in summary["workloads"]["fit"]
    ]
    workloads = ["fit", "selective", "heldout", "scan"]
    x = np.arange(len(methods))
    width = 0.18
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for i, workload in enumerate(workloads):
        vals = [summary["workloads"][workload][m]["ratio_to_generic_median"] for m in methods]
        ax.bar(x + (i - 1.5) * width, vals, width, label=workload)
    ax.axhline(1.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=25, ha="right")
    ax.set_ylabel("Objective / generic oracle")
    ax.set_title("Block-simulator objective relative to generic 2D tiling")
    ax.legend(frameon=False, ncol=4)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
