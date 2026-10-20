"""Query sets for fit, scan, selective, held-out, and short-window workloads."""

from __future__ import annotations

import numpy as np

from pedpack.grouping import nuclear_groups
from pedpack.pedigree import trios as all_trios
from pedpack.types import Cohort, Query


def _clip_window(n_variants: int, start: int, span: int) -> tuple[int, int]:
    span = max(1, min(span, n_variants))
    start = max(0, min(start, n_variants - span))
    return start, start + span


def _trios_inside(
    trios: list[tuple[int, int, int]],
    samples: set[int],
) -> tuple[tuple[int, int, int], ...]:
    return tuple(t for t in trios if t[0] in samples and t[1] in samples and t[2] in samples)


def make_workloads(
    cohort: Cohort,
    rng: np.random.Generator,
    n_fit_scan: int = 4,
    n_fit_selective: int = 12,
    n_eval_scan: int = 4,
    n_eval_selective: int = 8,
    n_heldout: int = 8,
    n_windows: int = 6,
    holdout_fraction: float = 0.30,
) -> dict[str, list[Query]]:
    """Build disjoint query lists.

    Nuclear families are split into a fit pool and a held-out pool. Scan
    queries read every sample. Selective queries read one nuclear family.
    Window queries are the short spans used for a compound-het style access
    pattern; they are not a caller.
    """

    families = nuclear_groups(cohort)
    order = rng.permutation(len(families))
    n_hold = max(1, int(round(len(families) * holdout_fraction)))
    if len(families) >= 3:
        n_hold = min(n_hold, len(families) - 2)
    hold_idx = set(order[:n_hold].tolist())
    train_idx = [int(i) for i in order.tolist() if int(i) not in hold_idx]
    hold_ids = [int(i) for i in order.tolist() if int(i) in hold_idx]
    if not train_idx:
        train_idx = hold_ids[:1]
    trio_list = all_trios(cohort.father, cohort.mother)
    n = cohort.n_samples
    m = cohort.n_variants
    everyone = tuple(range(n))
    all_trio_t = tuple(trio_list)

    def selective(name: str, workload: str, family_indexes: list[int], span_low: float, span_high: float, count: int) -> list[Query]:
        queries: list[Query] = []
        if not family_indexes or count <= 0:
            return queries
        for i in range(count):
            members = families[int(family_indexes[i % len(family_indexes)])]
            sample_set = set(members)
            frac = float(rng.uniform(span_low, span_high))
            span = max(1, int(round(frac * m)))
            start = int(rng.integers(0, max(1, m - span + 1)))
            v0, v1 = _clip_window(m, start, span)
            queries.append(
                Query(
                    name=f"{name}-{i}",
                    workload=workload,
                    samples=tuple(sorted(sample_set)),
                    v0=v0,
                    v1=v1,
                    trios=_trios_inside(trio_list, sample_set),
                )
            )
        return queries

    def scans(name: str, workload: str, span_low: float, span_high: float, count: int) -> list[Query]:
        queries: list[Query] = []
        for i in range(count):
            frac = float(rng.uniform(span_low, span_high))
            span = max(1, int(round(frac * m)))
            start = int(rng.integers(0, max(1, m - span + 1)))
            v0, v1 = _clip_window(m, start, span)
            queries.append(
                Query(
                    name=f"{name}-{i}",
                    workload=workload,
                    samples=everyone,
                    v0=v0,
                    v1=v1,
                    trios=all_trio_t,
                )
            )
        return queries

    def windows(count: int) -> list[Query]:
        queries: list[Query] = []
        if not hold_ids:
            return queries
        for i in range(count):
            members = families[hold_ids[i % len(hold_ids)]]
            sample_set = set(members)
            span = int(rng.integers(20, max(21, min(80, m))))
            start = int(rng.integers(0, max(1, m - span + 1)))
            v0, v1 = _clip_window(m, start, span)
            queries.append(
                Query(
                    name=f"window-{i}",
                    workload="window",
                    samples=tuple(sorted(sample_set)),
                    v0=v0,
                    v1=v1,
                    trios=_trios_inside(trio_list, sample_set),
                )
            )
        return queries

    fit = scans("fit-scan", "fit", 0.45, 1.0, n_fit_scan)
    fit += selective("fit-sel", "fit", train_idx, 0.05, 0.25, n_fit_selective)
    return {
        "fit": fit,
        "scan": scans("scan", "scan", 0.35, 1.0, n_eval_scan),
        "selective": selective("sel", "selective", train_idx, 0.08, 0.30, n_eval_selective),
        "heldout": selective("hold", "heldout", hold_ids, 0.01, 0.08, n_heldout // 2)
        + selective("hold-wide", "heldout", hold_ids, 0.40, 0.70, n_heldout - n_heldout // 2),
        "window": windows(n_windows),
    }
