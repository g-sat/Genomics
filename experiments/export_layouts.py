"""Write the layouts chosen on the positive regimes, plus the variant and generic ones.

Reads the already computed study JSON so the files match those choices.
Does not contact HDFS.

    python experiments/export_layouts.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pedpack.genotype import build_cohort  # noqa: E402
from pedpack.materialize import read_plan, write_layout  # noqa: E402
from pedpack.models import train_embeddings  # noqa: E402
from pedpack.pedigree import simulate_pedigree  # noqa: E402
from pedpack.study import _layouts, make_cluster  # noqa: E402
from pedpack.types import Query  # noqa: E402
from pedpack.workloads import make_workloads  # noqa: E402

SCENARIOS = ("select_05", "select_15", "mix_50")
METHODS = ("pedpack", "pedpack_no_pedigree_features", "generic_oracle")


def _queries_of(workloads: dict, names: tuple[str, ...]) -> list[Query]:
    out: list[Query] = []
    for name in names:
        out.extend(workloads.get(name, []))
    return out


def export_scenario(name: str) -> None:
    study = json.loads((ROOT / "results" / "stress" / name / "study.json").read_text(encoding="utf-8"))
    config = study["config"]
    cluster = make_cluster(config.get("cluster"))
    out_root = ROOT / "results" / "transfer" / name
    for seed_row in study["seeds"]:
        seed = int(seed_row["seed"])
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
        embeddings = train_embeddings(
            cohort,
            steps=int(config.get("encoder_steps", 80)),
            hidden=int(config.get("hidden", 32)),
            seed=seed,
        )
        layouts, _heuristic = _layouts(cohort, embeddings, cluster, config)
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
            n_mixed_scan=int(config.get("n_mixed_scan", 0)),
            n_mixed_family=int(config.get("n_mixed_family", 0)),
        )
        queries = _queries_of(workloads, ("family_selective", "heldout_family"))
        for method in METHODS:
            choice = seed_row["methods"][method]
            layout = layouts[int(choice["index"])]
            if layout.stats.label != choice["choice"]["label"]:
                raise RuntimeError(
                    f"{name} seed {seed} {method}: pool label {layout.stats.label} "
                    f"does not match stored choice {choice['choice']['label']}"
                )
            directory = out_root / f"seed_{seed}" / method
            manifest = write_layout(layout, cohort, directory)
            plans = []
            for query in queries:
                for item in read_plan(manifest, query):
                    item["workload"] = query.workload
                    plans.append(item)
            (directory / "reads.json").write_text(json.dumps(plans), encoding="utf-8")
            print(
                f"{name} seed {seed} {method}: {manifest['n_files']} files, "
                f"{manifest['primary_bytes']} bytes",
                flush=True,
            )


def main() -> None:
    names = sys.argv[1:] or list(SCENARIOS)
    for name in names:
        export_scenario(name)


if __name__ == "__main__":
    main()
