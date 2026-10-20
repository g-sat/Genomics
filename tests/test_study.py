"""One small study, checked for protocol invariants rather than a winner."""

from __future__ import annotations

import json
from pathlib import Path

from pedpack.study import run_study

ROOT = Path(__file__).resolve().parents[1]


def test_tiny_study_respects_the_measurement_budget():
    config = json.loads((ROOT / "configs" / "tiny.json").read_text(encoding="utf-8"))
    report = run_study(config, verbose=False)
    seed = report["seeds"][0]
    budget = config["budget"]
    assert seed["mendelian_errors"] == 0
    assert seed["methods"]["pedpack"]["layouts_measured"] <= budget
    assert seed["methods"]["analytical_budget"]["layouts_measured"] <= budget
    assert seed["methods"]["pedpack_generic_pool"]["layouts_measured"] <= budget
    assert seed["methods"]["random_budget"]["layouts_measured"] <= budget
    full = seed["methods"]["full_oracle"]["workloads"]["fit"]["objective"]
    pedpack = seed["methods"]["pedpack"]["workloads"]["fit"]["objective"]
    generic = seed["methods"]["generic_oracle"]["workloads"]["fit"]["objective"]
    assert full <= pedpack + 1e-6
    assert full <= generic + 1e-6
    assert seed["methods"]["heuristic_fixed"]["choice"]["trio_colocation"] == 1.0
    # The fit-tuned selector is scored on queries it was not allowed to fit.
    for workload in ("heldout", "window", "scan", "selective"):
        assert seed["methods"]["pedpack"]["workloads"][workload]["objective"] > 0
    assert seed["analysis"]["best_on_heldout"]["analysis_only"] is True
