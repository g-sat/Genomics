"""Selectivity of family-centric queries."""

from __future__ import annotations

import numpy as np

from pedpack.genotype import build_cohort
from pedpack.pedigree import simulate_pedigree
from pedpack.workloads import family_count, make_workloads


def test_point_fraction_rounds_half_up():
    rng = np.random.default_rng(0)
    assert family_count(10, 0.05, rng) == 1
    assert family_count(10, 0.15, rng) == 2
    assert family_count(10, 0.40, rng) == 4


def test_locked_draw_stays_inside_five_to_fifteen_percent():
    rng = np.random.default_rng(1)
    draws = [family_count(20, None, rng) for _ in range(30)]
    assert min(draws) >= 1
    assert max(draws) <= 3


def test_exact_fraction_changes_the_sample_set():
    rng = np.random.default_rng(3)
    pedigree = simulate_pedigree(10, 2, 2.2, 0.25, rng, min_people=40)
    cohort = build_cohort(pedigree, 80, rng)
    narrow = make_workloads(cohort, np.random.default_rng(4), family_fraction=0.05, n_family_selective=6, n_fit_family=2, n_fit_extended=1, n_extended=1, n_heldout_family=2, n_heldout_extended=1, n_fit_scan=1, n_eval_scan=1)
    wide = make_workloads(cohort, np.random.default_rng(4), family_fraction=0.40, n_family_selective=6, n_fit_family=2, n_fit_extended=1, n_extended=1, n_heldout_family=2, n_heldout_extended=1, n_fit_scan=1, n_eval_scan=1)
    narrow_n = np.mean([len(query.samples) for query in narrow["family_selective"]])
    wide_n = np.mean([len(query.samples) for query in wide["family_selective"]])
    assert wide_n > narrow_n
