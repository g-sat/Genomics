"""Re-fit the locked regime under selectivity, mixture, and pedigree changes.

The locked result in results/pedigree/ is not overwritten. Each scenario
starts from configs/pedigree_gap.json and changes one factor.

    python experiments/run_stress.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pedpack.report import write_report  # noqa: E402
from pedpack.study import run_study  # noqa: E402

SCENARIOS = (
    {
        "name": "select_05",
        "factor": "5% of families",
        "overlay": {"family_fraction": 0.05},
    },
    {
        "name": "select_15",
        "factor": "15% of families",
        "overlay": {"family_fraction": 0.15},
    },
    {
        "name": "select_40",
        "factor": "40% of families",
        "overlay": {"family_fraction": 0.40},
    },
    {
        "name": "mix_50",
        "factor": "fit mix 50% scan / 50% family",
        "overlay": {
            "fit_mixture": {"scan": 0.50, "family": 0.50, "extended": 0.0},
            "n_fit_extended": 0,
        },
    },
    {
        "name": "overlap",
        "factor": "3 generations, cross rate 0.6",
        "overlay": {"generations": 3, "cross_rate": 0.6, "min_people": 60},
    },
    {
        "name": "scale",
        "factor": "16 families, 3 generations, 3200 variants",
        "overlay": {
            "n_families": 16,
            "generations": 3,
            "cross_rate": 0.45,
            "min_people": 80,
            "n_variants": 3200,
            "widths": [200, 800, 1600],
            "heuristic_width": 800,
        },
    },
)


def _gap(report: dict) -> dict:
    rows = []
    for seed in report["seeds"]:
        ped = seed["methods"]["pedpack"]
        blind = seed["methods"]["pedpack_no_pedigree_features"]
        generic = seed["methods"]["generic_oracle"]
        variant = seed["methods"]["variant_oracle"]
        ped_cost = ped["workloads"]["fit"]["objective"]
        blind_cost = blind["workloads"]["fit"]["objective"]
        rows.append(
            {
                "seed": seed["seed"],
                "n_samples": seed["n_samples"],
                "n_trios": seed["n_trios"],
                "pedpack_family": ped["choice"]["family"],
                "blind_family": blind["choice"]["family"],
                "same_layout": ped["index"] == blind["index"],
                "pedpack": ped_cost,
                "blind": blind_cost,
                "variant": variant["workloads"]["fit"]["objective"],
                "generic": generic["workloads"]["fit"]["objective"],
                "blind_over_pedpack": blind_cost / ped_cost,
                "sample_fraction_family": seed["query_stats"]["family_selective"]["mean_sample_fraction"],
                "heldout_bytes_pedpack": ped["workloads"]["heldout_family"]["mean_bytes"],
                "heldout_bytes_blind": blind["workloads"]["heldout_family"]["mean_bytes"],
            }
        )
    return {"rows": rows, "gap_on_every_seed": all(not row["same_layout"] and row["blind_over_pedpack"] > 1.0 for row in rows)}


def _fmt(value: float) -> str:
    return f"{value:.3f}"


def _markdown(base: dict, blocks: list[dict]) -> str:
    lines = [
        "# Stress tests around the locked pedigree-gap regime",
        "",
        "Each scenario re-runs the search. Layouts are chosen on that scenario's",
        "fit mix, then scored on held-out families. Numbers are from the block",
        "simulator. They are not HDFS measurements.",
        "",
        f"Base config: seeds {base['seeds']}, horizon {base['horizon_queries']}, "
        f"mixture {base['fit_mixture']}, block {base['cluster']['block_size']}, "
        f"bytes per genotype {base['cluster']['bytes_per_genotype']}.",
        "",
        "| Scenario | Samples | PedPack | No-pedigree | Blind / PedPack | Same layout? | Held-out byte ratio | Gap |",
        "| --- | ---: | --- | --- | ---: | --- | ---: | --- |",
    ]
    for block in blocks:
        for row in block["gap"]["rows"]:
            byte_ratio = row["heldout_bytes_blind"] / row["heldout_bytes_pedpack"]
            lines.append(
                f"| {block['name']} | {row['n_samples']} | {row['pedpack_family']} | {row['blind_family']} | "
                f"{_fmt(row['blind_over_pedpack'])} | {row['same_layout']} | {_fmt(byte_ratio)} | "
                f"{'yes' if (not row['same_layout'] and row['blind_over_pedpack'] > 1) else 'no'} |"
            )
    lines.extend(["", "## What each scenario changes", ""])
    for block in blocks:
        held = "yes" if block["gap"]["gap_on_every_seed"] else "no"
        stats = block["report"]["seeds"][0]["query_stats"]["family_selective"]
        lines.append(
            f"- `{block['name']}`: {block['factor']}. "
            f"Family-selective queries touch {stats['mean_sample_fraction']:.1%} of samples "
            f"on seed {block['report']['seeds'][0]['seed']}. "
            f"Gap on every seed: {held}."
        )
    lines.extend(
        [
            "",
            "A gap means PedPack and the no-pedigree ranker chose different layouts and",
            "PedPack's measured fit objective is lower. The held-out byte ratio is",
            "no-pedigree bytes divided by PedPack bytes on families that were not used",
            "to choose the layout.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    base = json.loads((ROOT / "configs" / "pedigree_gap.json").read_text(encoding="utf-8"))
    out = ROOT / "results" / "stress"
    out.mkdir(parents=True, exist_ok=True)
    blocks = []
    for spec in SCENARIOS:
        config = dict(base)
        config.update(spec["overlay"])
        config["secondary_neural"] = False
        config["n_mixed_scan"] = 4
        config["n_mixed_family"] = 4
        print(f"scenario {spec['name']}", flush=True)
        report = run_study(config, verbose=True)
        write_report(report, out / spec["name"])
        gap = _gap(report)
        blocks.append({"name": spec["name"], "factor": spec["factor"], "report": report, "gap": gap})
        print(f"  gap_on_every_seed={gap['gap_on_every_seed']}", flush=True)
    text = _markdown(base, blocks)
    (out / "summary.md").write_text(text, encoding="utf-8")
    (out / "summary.json").write_text(
        json.dumps(
            [
                {
                    "name": block["name"],
                    "factor": block["factor"],
                    "gap_on_every_seed": block["gap"]["gap_on_every_seed"],
                    "rows": block["gap"]["rows"],
                }
                for block in blocks
            ],
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"wrote {out / 'summary.md'}")


if __name__ == "__main__":
    main()
