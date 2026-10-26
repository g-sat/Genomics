"""Co-located genotype and phenotype tiles.

Vectors are coordinates the caller supplies. This module does not read
clinical notes, an EHR, or a diagnostic platform. Similarity is cosine
similarity on those coordinates. Samples linked at or above a threshold form
one phenotype component.

A hyper-tile stores both payloads for the people in one component, so an
open pays for genotypes and phenotypes together. The split layout keeps
genotypes in pedigree groups and phenotypes in similarity components.
``choose_hyper_layout`` picks the cheaper of pedigree tiles, phenotype
tiles, and the split for a mix of Mendelian and differential queries.
"""

from __future__ import annotations

from pedpack.multimodal import group_bytes


def cosine(left: tuple[float, ...], right: tuple[float, ...]) -> float:
    if len(left) != len(right) or not left:
        raise ValueError("embeddings must share a positive dimension")
    dot = sum(a * b for a, b in zip(left, right, strict=True))
    left_norm = sum(a * a for a in left) ** 0.5
    right_norm = sum(b * b for b in right) ** 0.5
    if left_norm == 0.0 or right_norm == 0.0:
        raise ValueError("embedding must be non-zero")
    return dot / (left_norm * right_norm)


def phenotype_groups(
    vectors: list[tuple[float, ...]],
    threshold: float,
) -> list[list[int]]:
    """Connected components of the cosine graph."""

    if not 0.0 <= threshold <= 1.0:
        raise ValueError("cosine threshold must lie in [0, 1]")
    count = len(vectors)
    parent = list(range(count))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for left in range(count):
        for right in range(left + 1, count):
            if cosine(vectors[left], vectors[right]) >= threshold:
                union(left, right)
    grouped: dict[int, list[int]] = {}
    for index in range(count):
        grouped.setdefault(find(index), []).append(index)
    groups = [sorted(members) for members in grouped.values()]
    groups.sort(key=lambda members: members[0])
    return groups


def hyper_bytes(
    groups: list[list[int]],
    query: list[int],
    genotype_bytes: int,
    phenotype_bytes: int,
) -> int:
    """Opaque read of every opened tile, both payloads included."""

    if genotype_bytes < 1 or phenotype_bytes < 1:
        raise ValueError("both payloads must be positive")
    return group_bytes(groups, query, genotype_bytes + phenotype_bytes)


def split_bytes(
    pedigree: list[list[int]],
    phenotype: list[list[int]],
    query: list[int],
    genotype_bytes: int,
    phenotype_bytes: int,
    mode: str,
) -> int:
    """Separate files. ``joint`` opens both groupings for the same ids."""

    if mode == "mendelian":
        return group_bytes(pedigree, query, genotype_bytes)
    if mode == "differential":
        return group_bytes(phenotype, query, phenotype_bytes)
    if mode == "joint":
        return group_bytes(pedigree, query, genotype_bytes) + group_bytes(
            phenotype, query, phenotype_bytes
        )
    raise ValueError("mode must be mendelian, differential, or joint")


def mix_bytes(fraction: float, mendelian: int, differential: int) -> float:
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("Mendelian fraction must lie in [0, 1]")
    return fraction * mendelian + (1.0 - fraction) * differential


def choose_hyper_layout(
    fraction: float,
    options: dict[str, tuple[int, int]],
) -> str:
    """Pick the layout with the lowest mixed byte cost.

    Each option maps a name to ``(Mendelian bytes, differential bytes)``.
    """

    if not options:
        raise ValueError("at least one layout is required")
    scored = {
        name: mix_bytes(fraction, mendelian, differential)
        for name, (mendelian, differential) in options.items()
    }
    best = min(scored.values())
    winners = [name for name, score in scored.items() if score == best]
    if len(winners) == 1:
        return winners[0]
    return "tie"
