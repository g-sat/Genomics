"""Query sets that can make an all-samples layout the wrong choice.

A scan still reads every sample. A family-centric query reads 5–15% of the
nuclear families on one genomic window. An extended query reads a nuclear
family plus the parents of its members, so grandparents and a shared parent
generation are in the request even when they are not in that nuclear group.
Fit and held-out lists use different families and different windows.
"""

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


def _round_half_up(value: float) -> int:
    return int(np.floor(float(value) + 0.5))


def family_count(n_families: int, fraction: float | tuple[float, float] | None, rng: np.random.Generator) -> int:
    """How many nuclear families one selective query touches.

    ``None`` keeps the locked regime: a uniform draw from 5% to 15% of the
    families, at least one. A single fraction is an exact count (half rounds
    up, and the count is at least one). A pair is an inclusive range.
    """

    if n_families <= 1:
        return 1
    if fraction is None:
        low = max(1, int(np.ceil(0.05 * n_families)))
        high = max(low, int(np.ceil(0.15 * n_families)))
        return int(rng.integers(low, high + 1))
    if isinstance(fraction, (int, float)):
        return int(min(n_families, max(1, _round_half_up(float(fraction) * n_families))))
    low_f, high_f = float(fraction[0]), float(fraction[1])
    low = max(1, _round_half_up(low_f * n_families))
    high = max(low, _round_half_up(high_f * n_families))
    high = min(n_families, high)
    low = min(low, high)
    if low == high:
        return low
    return int(rng.integers(low, high + 1))


def _extended_members(cohort: Cohort, members: list[int] | tuple[int, ...]) -> list[int]:
    """Nuclear members plus their parents. Parents of founders are ignored."""

    extra = set(int(i) for i in members)
    for person in list(extra):
        for parent in (int(cohort.father[person]), int(cohort.mother[person])):
            if parent >= 0:
                extra.add(parent)
    return sorted(extra)


def make_workloads(
    cohort: Cohort,
    rng: np.random.Generator,
    n_fit_scan: int = 4,
    n_fit_family: int = 8,
    n_fit_extended: int = 6,
    n_eval_scan: int = 4,
    n_family_selective: int = 8,
    n_extended: int = 6,
    n_heldout_family: int = 6,
    n_heldout_extended: int = 4,
    holdout_fraction: float = 0.30,
    family_span: tuple[float, float] = (0.20, 0.40),
    extended_span: tuple[float, float] = (0.12, 0.28),
    scan_span: tuple[float, float] = (0.55, 1.0),
    family_fraction: float | tuple[float, float] | None = None,
    n_mixed_scan: int = 0,
    n_mixed_family: int = 0,
) -> dict[str, list[Query]]:
    """Build the fit mixture and the evaluation lists.

    Held-out queries are seeded on families that do not appear as the seed of
    a fit query. An extended held-out query may still name a parent who also
    sits in a training family; that person is a relative, not a second seed.
    """

    families = nuclear_groups(cohort)
    order = rng.permutation(len(families))
    n_hold = max(1, int(round(len(families) * holdout_fraction)))
    if len(families) >= 3:
        n_hold = min(n_hold, len(families) - 2)
    hold_idx = set(int(i) for i in order[:n_hold].tolist())
    train_idx = [int(i) for i in order.tolist() if int(i) not in hold_idx]
    hold_ids = [int(i) for i in order.tolist() if int(i) in hold_idx]
    if not train_idx:
        train_idx = hold_ids[:1]
    trio_list = all_trios(cohort.father, cohort.mother)
    n = cohort.n_samples
    m = cohort.n_variants
    everyone = tuple(range(n))
    all_trio_t = tuple(trio_list)

    def window(span_low: float, span_high: float) -> tuple[int, int]:
        frac = float(rng.uniform(span_low, span_high))
        span = max(1, int(round(frac * m)))
        start = int(rng.integers(0, max(1, m - span + 1)))
        return _clip_window(m, start, span)

    def scans(name: str, workload: str, count: int) -> list[Query]:
        queries: list[Query] = []
        for i in range(count):
            v0, v1 = window(*scan_span)
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

    def family_queries(name: str, workload: str, pool: list[int], count: int) -> list[Query]:
        queries: list[Query] = []
        if not pool or count <= 0:
            return queries
        for i in range(count):
            k = min(family_count(len(families), family_fraction, rng), len(pool))
            chosen = rng.choice(pool, size=k, replace=False)
            members: set[int] = set()
            for index in np.atleast_1d(chosen).tolist():
                members.update(int(person) for person in families[int(index)])
            v0, v1 = window(*family_span)
            queries.append(
                Query(
                    name=f"{name}-{i}",
                    workload=workload,
                    samples=tuple(sorted(members)),
                    v0=v0,
                    v1=v1,
                    trios=_trios_inside(trio_list, members),
                )
            )
        return queries

    def extended_queries(name: str, workload: str, pool: list[int], count: int) -> list[Query]:
        queries: list[Query] = []
        if not pool or count <= 0:
            return queries
        growing = [i for i in pool if len(_extended_members(cohort, families[i])) > len(families[i])]
        source = growing or pool
        for i in range(count):
            seed = int(source[int(rng.integers(0, len(source)))])
            members = set(_extended_members(cohort, families[seed]))
            # A second overlapping family, when one exists, shares a parent generation.
            others = [j for j in pool if j != seed]
            if others and float(rng.random()) < 0.5:
                extra = int(others[int(rng.integers(0, len(others)))])
                members.update(_extended_members(cohort, families[extra]))
            v0, v1 = window(*extended_span)
            queries.append(
                Query(
                    name=f"{name}-{i}",
                    workload=workload,
                    samples=tuple(sorted(members)),
                    v0=v0,
                    v1=v1,
                    trios=_trios_inside(trio_list, members),
                )
            )
        return queries

    fit = scans("fit-scan", "fit", n_fit_scan)
    fit += family_queries("fit-family", "fit", train_idx, n_fit_family)
    fit += extended_queries("fit-ext", "fit", train_idx, n_fit_extended)
    out = {
        "fit": fit,
        "scan": scans("scan", "scan", n_eval_scan),
        "family_selective": family_queries("fam", "family_selective", train_idx, n_family_selective),
        "extended": extended_queries("ext", "extended", train_idx, n_extended),
        "heldout_family": family_queries("hold-fam", "heldout_family", hold_ids, n_heldout_family),
        "heldout_extended": extended_queries("hold-ext", "heldout_extended", hold_ids, n_heldout_extended),
    }
    # Drawn last so the locked fit and held-out lists keep their random stream.
    mixed = scans("mixed-scan", "mixed", n_mixed_scan)
    mixed += family_queries("mixed-family", "mixed", train_idx, n_mixed_family)
    if mixed:
        out["mixed"] = mixed
    return out
