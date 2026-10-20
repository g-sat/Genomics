"""Budgeted layout selection.

Oracles are not implemented here. They scan every candidate in a family and
are reported separately so a reader can tell a search with a measurement
budget from an exhaustive table.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pedpack.cost import analytical_query, objective
from pedpack.features import feature_matrix, stack_targets
from pedpack.models import Surrogate, train_surrogate
from pedpack.types import Cluster, Cohort, Measurement, PreparedLayout, Query


@dataclass
class Selection:
    """Result of a budgeted search. ``measured`` counts layouts, not queries."""

    index: int
    measured: int
    pilot: int
    shortlist: int
    surrogate: Surrogate | None
    predicted_objective: np.ndarray | None
    notes: str
    measured_indexes: tuple[int, ...] = ()


def _weighted_seconds(seconds: list[float] | np.ndarray, query_weights: np.ndarray | None) -> float:
    values = np.asarray(seconds, dtype=np.float64)
    if query_weights is None:
        return float(values.sum())
    weights = np.asarray(query_weights, dtype=np.float64)
    if weights.shape[0] != values.shape[0]:
        raise ValueError("query weights do not match the query list")
    return float(np.dot(weights, values))


def measured_objective(
    layout: PreparedLayout,
    measurements: list[Measurement],
    weight: float,
    query_weights: np.ndarray | None = None,
) -> float:
    if any(m.uncovered for m in measurements):
        return 1.0e12
    seconds = _weighted_seconds([m.seconds for m in measurements], query_weights)
    return objective(layout.stats, seconds, weight)


def analytical_objective(
    layout: PreparedLayout,
    queries: list[Query],
    cluster: Cluster,
    weight: float,
    query_weights: np.ndarray | None = None,
    pedigree_aware: bool = True,
) -> float:
    seconds: list[float] = []
    for query in queries:
        result = analytical_query(layout, query, cluster, pedigree_aware=pedigree_aware)
        if result.uncovered:
            return 1.0e12
        seconds.append(result.seconds)
    return objective(layout.stats, _weighted_seconds(seconds, query_weights), weight)


def pilot_count(n_layouts: int, budget: int, pilot_fraction: float) -> int:
    """How many layouts the neural selector measures before it trusts the network."""

    take = min(budget, n_layouts)
    n_pilot = int(round(take * pilot_fraction))
    if take >= 6:
        return min(max(4, n_pilot), take - 1)
    return take


def stratified_pilot(layouts: list[PreparedLayout], count: int, rng: np.random.Generator) -> list[int]:
    """Draw a pilot that touches variant, generic, and pedigree layouts."""

    buckets: dict[str, list[int]] = {"variant": [], "generic": [], "pedigree": []}
    for index, layout in enumerate(layouts):
        family = layout.stats.family
        key = family if family in {"variant", "generic"} else "pedigree"
        buckets[key].append(index)
    for values in buckets.values():
        rng.shuffle(values)
    if count >= len(layouts):
        return list(range(len(layouts)))
    chosen: list[int] = []
    while len(chosen) < count and any(buckets.values()):
        for key in ("variant", "generic", "pedigree"):
            if buckets[key] and len(chosen) < count:
                chosen.append(buckets[key].pop())
    if len(chosen) < count:
        rest = [i for i in range(len(layouts)) if i not in set(chosen)]
        rng.shuffle(rest)
        chosen.extend(rest[: count - len(chosen)])
    return chosen


def select_random(
    layouts: list[PreparedLayout],
    queries: list[Query],
    cluster: Cluster,
    weight: float,
    budget: int,
    rng: np.random.Generator,
    measure,
    query_weights: np.ndarray | None = None,
) -> Selection:
    del queries, cluster
    take = min(budget, len(layouts))
    indexes = rng.choice(len(layouts), size=take, replace=False).tolist()
    best_i = indexes[0]
    best_cost = float("inf")
    for index in indexes:
        cost = measured_objective(layouts[index], measure(index), weight, query_weights)
        if cost < best_cost:
            best_cost = cost
            best_i = index
    return Selection(best_i, take, 0, take, None, None, "random pilot")


def select_analytical(
    layouts: list[PreparedLayout],
    queries: list[Query],
    cluster: Cluster,
    weight: float,
    budget: int,
    measure,
    query_weights: np.ndarray | None = None,
    pedigree_aware: bool = True,
) -> Selection:
    scores = np.asarray(
        [
            analytical_objective(
                layout, queries, cluster, weight, query_weights, pedigree_aware=pedigree_aware
            )
            for layout in layouts
        ],
        dtype=np.float64,
    )
    # mergesort is stable, so equal scores keep candidate order instead of hash order
    order = np.argsort(scores, kind="mergesort")
    take = min(budget, len(layouts))
    shortlist = [int(i) for i in order[:take].tolist()]
    best_i = shortlist[0]
    best_cost = float("inf")
    for index in shortlist:
        cost = measured_objective(layouts[index], measure(index), weight, query_weights)
        if cost < best_cost:
            best_cost = cost
            best_i = index
    note = "analytical top-B, pedigree-aware" if pedigree_aware else "analytical top-B, pedigree terms removed"
    return Selection(best_i, take, 0, take, None, scores, note, tuple(shortlist))


def select_neural(
    layouts: list[PreparedLayout],
    queries: list[Query],
    cohort: Cohort,
    cluster: Cluster,
    embeddings: np.ndarray | None,
    weight: float,
    budget: int,
    pilot_fraction: float,
    rng: np.random.Generator,
    measure,
    hide_pedigree: bool = False,
    hidden: int = 32,
    epochs: int = 80,
    seed: int = 0,
    pilot_ids: list[int] | None = None,
    query_weights: np.ndarray | None = None,
) -> Selection:
    """Fit on a stratified pilot, remeasure the predicted top of the remainder."""

    take = min(budget, len(layouts))
    n_pilot = pilot_count(len(layouts), budget, pilot_fraction)
    pilot = list(pilot_ids) if pilot_ids is not None else stratified_pilot(layouts, n_pilot, rng)
    pilot_set = set(pilot)
    train_measures: list[Measurement] = []
    train_layouts: list[PreparedLayout] = []
    for index in pilot:
        train_layouts.append(layouts[index])
        train_measures.extend(measure(index))
    # train_measures is grouped by layout; that matches feature_matrix, which is layout-major.
    features = feature_matrix(
        train_layouts,
        queries,
        cohort,
        cluster,
        embeddings,
        hide_pedigree=hide_pedigree,
    )
    targets = stack_targets(train_measures)
    if features.shape[0] != targets.shape[0]:
        raise RuntimeError("training rows are not aligned with measurements")
    surrogate = train_surrogate(features, targets, hidden=hidden, epochs=epochs, seed=seed)
    scores = np.full(len(layouts), np.inf, dtype=np.float64)
    for index in pilot:
        scores[index] = measured_objective(layouts[index], measure(index), weight, query_weights)
    # The second measure() call hits the cache. Count it once, below.
    complement = [i for i in range(len(layouts)) if i not in pilot_set]
    shortlist: list[int] = []
    if complement:
        pred_matrix = predict_query_seconds(
            surrogate,
            [layouts[i] for i in complement],
            queries,
            cohort,
            cluster,
            embeddings,
            hide_pedigree,
        )
        pred_seconds = _combine_predictions(pred_matrix, query_weights)
        pred_obj = np.asarray(
            [
                objective(layouts[index].stats, float(pred_seconds[pos]), weight)
                for pos, index in enumerate(complement)
            ],
            dtype=np.float64,
        )
        for pos, index in enumerate(complement):
            scores[index] = float(pred_obj[pos])
        room = take - len(pilot)
        if room > 0:
            order = np.argsort(pred_obj)
            shortlist = [complement[int(pos)] for pos in order[:room].tolist()]
    best_i = pilot[0]
    best_cost = float("inf")
    considered = list(pilot) + shortlist
    seen: set[int] = set()
    for index in considered:
        if index in seen:
            continue
        seen.add(index)
        cost = measured_objective(layouts[index], measure(index), weight, query_weights)
        if cost < best_cost:
            best_cost = cost
            best_i = index
    note = "neural shortlist"
    if hide_pedigree:
        note = "neural shortlist, pedigree features hidden"
    return Selection(
        index=best_i,
        measured=len(seen),
        pilot=len(pilot),
        shortlist=len(shortlist),
        surrogate=surrogate,
        predicted_objective=scores,
        notes=note,
    )


def _combine_predictions(pred_matrix: np.ndarray, query_weights: np.ndarray | None) -> np.ndarray:
    if query_weights is None:
        return pred_matrix.sum(axis=1)
    weights = np.asarray(query_weights, dtype=np.float64)
    return pred_matrix @ weights


def predict_query_seconds(
    surrogate: Surrogate,
    layouts: list[PreparedLayout],
    queries: list[Query],
    cohort: Cohort,
    cluster: Cluster,
    embeddings: np.ndarray | None,
    hide_pedigree: bool,
) -> np.ndarray:
    """Predicted seconds with shape ``(n_layouts, n_queries)``."""

    features = feature_matrix(
        layouts, queries, cohort, cluster, embeddings, hide_pedigree=hide_pedigree
    )
    decoded = surrogate.predict(features)
    return decoded[:, 3].reshape(len(layouts), len(queries))
