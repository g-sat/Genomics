"""Choose a pedigree cut or a phenotype cut for one workload mix.

Phenotype groups are an assignment supplied by the caller. They are not
clinical notes. A hereditary query names one pedigree group. A phenotype
query names one phenotype group. The store reads every sample in each
opened group.

Let ``h`` be the fraction of queries that are hereditary. Bytes per query are

    C_ped = h * H_ped + (1 - h) * P_ped
    C_phe = h * H_phe + (1 - h) * P_phe

Pedigree is cheaper when ``h`` is above the threshold below, provided
pedigree is the better hereditary layout and phenotype is the better
phenotype layout.
"""

from __future__ import annotations


def group_bytes(groups: list[list[int]], query: list[int], sample_bytes: int) -> int:
    """Bytes to read every group that contains one queried sample."""

    if sample_bytes < 1:
        raise ValueError("sample payload must be positive")
    wanted = set(query)
    if not wanted:
        return 0
    total = 0
    covered = set()
    for group in groups:
        members = set(group)
        if members & wanted:
            total += len(group) * sample_bytes
            covered |= members & wanted
    if covered != wanted:
        raise ValueError("query samples are missing from the grouping")
    return total


def pedigree_threshold(h_ped: int, h_phe: int, p_ped: int, p_phe: int) -> float | None:
    """Hereditary fraction above which the pedigree cut is cheaper.

    ``None`` means one layout is cheaper at every fraction.
    """

    # h * h_ped + (1-h) * p_ped < h * h_phe + (1-h) * p_phe
    # h * (h_ped - p_ped - h_phe + p_phe) < p_phe - p_ped
    left = h_ped - p_ped - h_phe + p_phe
    right = p_phe - p_ped
    if left == 0:
        return None
    return right / left


def choose_cut(h: float, h_ped: int, h_phe: int, p_ped: int, p_phe: int) -> str:
    """Return ``pedigree`` or ``phenotype`` for a hereditary fraction ``h``."""

    if not 0.0 <= h <= 1.0:
        raise ValueError("hereditary fraction must lie in [0, 1]")
    pedigree = h * h_ped + (1.0 - h) * p_ped
    phenotype = h * h_phe + (1.0 - h) * p_phe
    if pedigree < phenotype:
        return "pedigree"
    if phenotype < pedigree:
        return "phenotype"
    return "tie"
