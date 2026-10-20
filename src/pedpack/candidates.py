"""Candidate layouts for the baselines and for PedPack's search space."""

from __future__ import annotations

import numpy as np

from pedpack.grouping import agglomerate, generic_chunks, nuclear_groups, three_generation_groups
from pedpack.types import Candidate, Cohort


def _as_groups(groups: list[list[int]]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(x) for x in sorted(group)) for group in groups)


def build_pool(
    cohort: Cohort,
    embeddings: np.ndarray | None,
    widths: list[int],
    generic_sizes: list[int],
    agg_sizes: list[int],
    min_kinship: list[float],
    lams: list[float],
    orders: list[str],
    packs: list[str],
    heuristic_width: int,
    heuristic_pack: str,
    heuristic_order: str,
    include_three_gen: bool = True,
) -> tuple[list[Candidate], int]:
    """Return the candidate list and the index of the fixed pedigree heuristic.

    Duplicate geometries are dropped. The heuristic index points at the fixed
    nuclear layout even when an identical candidate was already added under
    the ``family`` tag.
    """

    n = cohort.n_samples
    m = cohort.n_variants
    widths = sorted({max(1, min(int(w), m)) for w in widths})
    pack_list = list(packs)
    order_list = list(orders)
    pool: list[Candidate] = []
    seen: dict[tuple, int] = {}

    def add(family: str, groups: list[list[int]], width: int, pack: str, order: str) -> int:
        grouped = _as_groups(groups)
        key = (grouped, int(width), pack, order)
        if key in seen:
            return seen[key]
        label = f"{family}|g={len(grouped)}|w={width}|{pack}|{order}|id={len(pool)}"
        candidate = Candidate(
            label=label,
            family=family,
            groups=grouped,
            width=int(width),
            pack=pack,
            order=order,
        )
        seen[key] = len(pool)
        pool.append(candidate)
        return seen[key]

    def expand(family: str, groups: list[list[int]]) -> None:
        for width in widths:
            for pack in pack_list:
                for order in order_list:
                    add(family, groups, width, pack, order)

    expand("variant", [list(range(n))])
    for size in generic_sizes:
        expand("generic", generic_chunks(n, int(size)))

    nuclear = nuclear_groups(cohort)
    expand("family", nuclear)
    if include_three_gen:
        expand("pedigree", three_generation_groups(cohort, nuclear))
    for max_size in agg_sizes:
        for threshold in min_kinship:
            for lam in lams:
                if lam != 0.0 and embeddings is None:
                    continue
                groups = agglomerate(
                    cohort.kinship,
                    max_size=int(max_size),
                    min_kinship=float(threshold),
                    embeddings=embeddings,
                    lam=float(lam),
                )
                expand("pedigree", groups)

    h_width = max(1, min(int(heuristic_width), m))
    heuristic_index = add("heuristic", nuclear, h_width, heuristic_pack, heuristic_order)
    return pool, heuristic_index
