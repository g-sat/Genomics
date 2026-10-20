"""Sample groupings: generic chunks, nuclear families, and kinship agglomeration."""

from __future__ import annotations

import numpy as np

from pedpack.types import Cohort


def nuclear_groups(cohort: Cohort) -> list[list[int]]:
    """Parents plus their children. Childless people form singleton groups."""

    buckets: dict[tuple[int, int], set[int]] = {}
    n = cohort.n_samples
    for child in range(n):
        p = int(cohort.father[child])
        q = int(cohort.mother[child])
        if p < 0 or q < 0:
            continue
        key = (p, q) if p < q else (q, p)
        buckets.setdefault(key, set()).update((child, p, q))
    covered: set[int] = set()
    groups = []
    for members in buckets.values():
        groups.append(sorted(members))
        covered.update(members)
    for i in range(n):
        if i not in covered:
            groups.append([i])
    return groups


def three_generation_groups(cohort: Cohort, nuclear: list[list[int]] | None = None) -> list[list[int]]:
    """Nuclear families expanded by the parents of each member.

    A person can sit in more than one group. Storage pays for the duplicate;
    a query planner may read only one copy.
    """

    if nuclear is None:
        nuclear = nuclear_groups(cohort)
    uniq: dict[tuple[int, ...], None] = {}
    for group in nuclear:
        extra = set(group)
        for person in list(extra):
            for parent in (int(cohort.father[person]), int(cohort.mother[person])):
                if parent >= 0:
                    extra.add(parent)
        uniq[tuple(sorted(extra))] = None
    return [list(key) for key in uniq]


def generic_chunks(n_samples: int, size: int) -> list[list[int]]:
    """Index-contiguous chunks. Column order is scrambled, so this ignores kinship."""

    if size < 1:
        raise ValueError("chunk size must be positive")
    return [list(range(start, min(start + size, n_samples))) for start in range(0, n_samples, size)]


def _cross_kinship(kinship: np.ndarray, left: list[int], right: list[int]) -> float:
    block = kinship[np.ix_(left, right)]
    return float(block.mean()) if block.size else 0.0


def _mean_vector(embeddings: np.ndarray, members: list[int]) -> np.ndarray:
    return embeddings[members].mean(axis=0)


def agglomerate(
    kinship: np.ndarray,
    max_size: int,
    min_kinship: float,
    embeddings: np.ndarray | None = None,
    lam: float = 0.0,
) -> list[list[int]]:
    """Greedy merge of highest-scoring pairs that respect ``max_size``.

    The score is mean cross-kinship, plus ``lam`` times the cosine of the
    cluster embedding means when embeddings are provided. Ties break toward
    the earlier cluster index, so the result is deterministic.
    """

    n = int(kinship.shape[0])
    if max_size < 1:
        raise ValueError("max_size must be positive")
    if embeddings is not None and lam != 0.0:
        return _agglomerate_with_embeddings(kinship, max_size, min_kinship, embeddings, lam)
    return _agglomerate_kinship(kinship, max_size, min_kinship)


def _agglomerate_kinship(
    kinship: np.ndarray,
    max_size: int,
    min_kinship: float,
) -> list[list[int]]:
    """Average-linkage merges. The score update is the size-weighted mean.

    On dyadic kinship values this matches recomputing the block mean. The
    first maximum in row-major order of the upper triangle is the same pair
    the nested scan used to keep: smaller index first, and a strict improvement
    over ``min_kinship``.
    """

    n = int(kinship.shape[0])
    members: list[list[int]] = [[i] for i in range(n)]
    sizes = np.ones(n, dtype=np.int32)
    alive = np.ones(n, dtype=bool)
    score = np.array(kinship, dtype=np.float64, copy=True)
    np.fill_diagonal(score, -np.inf)

    while int(alive.sum()) > 1:
        masked = score.copy()
        too_big = (sizes[:, None] + sizes[None, :]) > max_size
        masked[too_big] = -np.inf
        masked[~alive, :] = -np.inf
        masked[:, ~alive] = -np.inf
        masked[np.tril_indices(n)] = -np.inf
        flat = int(np.argmax(masked))
        a, b = divmod(flat, n)
        best = float(masked[a, b])
        if not np.isfinite(best) or best <= min_kinship:
            break
        others = np.flatnonzero(alive)
        others = others[(others != a) & (others != b)]
        na = float(sizes[a])
        nb = float(sizes[b])
        updated = (na * score[a, others] + nb * score[b, others]) / (na + nb)
        score[a, others] = updated
        score[others, a] = updated
        members[a].extend(members[b])
        sizes[a] = int(na + nb)
        alive[b] = False
    return [sorted(members[i]) for i in range(n) if alive[i]]


def _agglomerate_with_embeddings(
    kinship: np.ndarray,
    max_size: int,
    min_kinship: float,
    embeddings: np.ndarray,
    lam: float,
) -> list[list[int]]:
    n = int(kinship.shape[0])
    clusters: list[list[int]] = [[i] for i in range(n)]
    means: list[np.ndarray] = [embeddings[i].copy() for i in range(n)]

    while len(clusters) > 1:
        best: tuple[int, int] | None = None
        best_score = min_kinship
        for a in range(len(clusters)):
            for b in range(a + 1, len(clusters)):
                if len(clusters[a]) + len(clusters[b]) > max_size:
                    continue
                score = _cross_kinship(kinship, clusters[a], clusters[b])
                score += lam * _cosine(means[a], means[b])
                if score > best_score:
                    best_score = score
                    best = (a, b)
        if best is None:
            break
        a, b = best
        clusters[a] = clusters[a] + clusters[b]
        wa = float(len(clusters[a]) - len(clusters[b]))
        wb = float(len(clusters[b]))
        means[a] = (means[a] * wa + means[b] * wb) / (wa + wb)
        del clusters[b]
        del means[b]
    return [sorted(c) for c in clusters]


def _cosine(left: np.ndarray, right: np.ndarray) -> float:
    denom = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denom < 1e-12:
        return 0.0
    return float(np.dot(left, right) / denom)


def parent_child_edges(cohort: Cohort) -> list[tuple[int, int]]:
    edges: list[tuple[int, int]] = []
    for child in range(cohort.n_samples):
        for parent in (int(cohort.father[child]), int(cohort.mother[child])):
            if parent >= 0:
                edges.append((parent, child))
    return edges


def grouping_stats(groups: list[list[int]] | tuple[tuple[int, ...], ...], cohort: Cohort) -> dict[str, float]:
    """Scalar descriptions used by both the analytical model and the ranker."""

    sets = [set(g) for g in groups]
    n = cohort.n_samples
    sizes = np.asarray([len(g) for g in groups], dtype=np.float64)
    duplication = float(sizes.sum() / n) if n else 0.0
    mean_size = float(sizes.mean()) if sizes.size else 0.0
    cv = float(sizes.std() / mean_size) if mean_size > 0 else 0.0

    trio_hit = 0
    trio_total = 0
    for child in range(n):
        p = int(cohort.father[child])
        q = int(cohort.mother[child])
        if p < 0 or q < 0:
            continue
        trio_total += 1
        if any(p in s and q in s and child in s for s in sets):
            trio_hit += 1
    trio_colocation = (trio_hit / trio_total) if trio_total else 1.0

    cut = 0
    edges = parent_child_edges(cohort)
    for parent, child in edges:
        if not any(parent in s and child in s for s in sets):
            cut += 1
    pc_cut = (cut / len(edges)) if edges else 0.0

    kin_num = 0.0
    kin_den = 0.0
    for group in groups:
        if len(group) < 2:
            continue
        block = cohort.kinship[np.ix_(group, group)]
        mask = ~np.eye(len(group), dtype=bool)
        kin_num += float(block[mask].mean()) * len(group)
        kin_den += len(group)
    within = (kin_num / kin_den) if kin_den else 0.0
    return {
        "n_groups": float(len(groups)),
        "mean_group_size": mean_size,
        "group_size_cv": cv,
        "trio_colocation": trio_colocation,
        "pc_cut_rate": pc_cut,
        "mean_within_kinship": within,
        "duplication": duplication,
    }
