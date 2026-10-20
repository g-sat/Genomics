"""Pedigree simulation, PED parsing, and kinship.

Kinship follows the standard recursive definition (Karigl / Lange): parents
have lower indices than their children, founders have self-kinship 1/2, and

    phi(i, i) = 1/2 * (1 + phi(p, q))
    phi(i, j) = 1/2 * (phi(j, p) + phi(j, q))

for j < i, with a missing parent contributing nothing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Pedigree:
    """Individuals in topological order. Parent pointers are indices or -1."""

    ids: np.ndarray
    father: np.ndarray
    mother: np.ndarray
    sex: np.ndarray
    family_id: np.ndarray

    @property
    def n(self) -> int:
        return int(self.ids.shape[0])

    @property
    def founder(self) -> np.ndarray:
        return (self.father < 0) & (self.mother < 0)


def kinship_matrix(father: np.ndarray, mother: np.ndarray) -> np.ndarray:
    """Return the symmetric kinship matrix for a topologically ordered pedigree."""

    father = np.asarray(father, dtype=np.int32)
    mother = np.asarray(mother, dtype=np.int32)
    n = int(father.shape[0])
    if mother.shape[0] != n:
        raise ValueError("father and mother vectors differ in length")
    phi = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        p = int(father[i])
        q = int(mother[i])
        if p < 0 and q < 0:
            phi[i, i] = 0.5
            continue
        if p >= i or q >= i:
            raise ValueError(
                f"individual {i} is not after its parents; sort the pedigree first"
            )
        related = 0.0 if p < 0 or q < 0 else float(phi[p, q])
        phi[i, i] = 0.5 * (1.0 + related)
        for j in range(i):
            s = 0.0
            if p >= 0:
                s += float(phi[j, p])
            if q >= 0:
                s += float(phi[j, q])
            phi[i, j] = phi[j, i] = 0.5 * s
    return phi


def trios(father: np.ndarray, mother: np.ndarray) -> list[tuple[int, int, int]]:
    """Each non-founder is one trio ``(father, mother, child)``."""

    out: list[tuple[int, int, int]] = []
    for child in range(int(father.shape[0])):
        p = int(father[child])
        q = int(mother[child])
        if p >= 0 and q >= 0:
            out.append((p, q, child))
    return out


def topological_order(
    ids: list[str],
    father_id: list[str | None],
    mother_id: list[str | None],
) -> np.ndarray:
    """Return positions into ``ids`` such that every parent precedes its child."""

    index = {ids[i]: i for i in range(len(ids))}
    n = len(ids)
    children: list[list[int]] = [[] for _ in range(n)]
    indeg = np.zeros(n, dtype=np.int32)
    for i in range(n):
        for parent in (father_id[i], mother_id[i]):
            if parent is None:
                continue
            if parent not in index:
                continue
            p = index[parent]
            children[p].append(i)
            indeg[i] += 1
    queue = [i for i in range(n) if indeg[i] == 0]
    order: list[int] = []
    while queue:
        node = queue.pop(0)
        order.append(node)
        for child in children[node]:
            indeg[child] -= 1
            if indeg[child] == 0:
                queue.append(child)
    if len(order) != n:
        raise ValueError("pedigree contains a cycle")
    return np.asarray(order, dtype=np.int32)


def parse_ped(path: str | Path) -> Pedigree:
    """Read a PLINK-style PED file (family, id, father, mother, sex, phenotype).

    Phenotype is ignored. Parent id ``0`` means a founder. A parent id that is
    absent from the file is treated as unknown, which is what a trimmed
    pedigree extract usually means.
    """

    ids: list[str] = []
    fams: list[str] = []
    father_ids: list[str | None] = []
    mother_ids: list[str | None] = []
    sexes: list[int] = []
    text = Path(path).read_text(encoding="utf-8")
    for line_no, line in enumerate(text.splitlines(), start=1):
        raw = line.strip()
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split()
        if len(parts) < 5:
            raise ValueError(f"PED line {line_no} has {len(parts)} fields, need at least 5")
        fam, iid, fid, mid, sex = parts[:5]
        ids.append(iid)
        fams.append(fam)
        father_ids.append(None if fid in {"0", "-9", ""} else fid)
        mother_ids.append(None if mid in {"0", "-9", ""} else mid)
        try:
            sexes.append(int(float(sex)))
        except ValueError as exc:
            raise ValueError(f"PED line {line_no} has a non-numeric sex") from exc
    if not ids:
        raise ValueError(f"no individuals in {path}")
    if len(set(ids)) != len(ids):
        raise ValueError("PED individual ids are not unique")
    known = set(ids)
    father_ids = [p if p in known else None for p in father_ids]
    mother_ids = [p if p in known else None for p in mother_ids]
    order = topological_order(ids, father_ids, mother_ids)
    id_of = {ids[i]: i for i in range(len(ids))}
    new_ids = np.asarray([ids[i] for i in order], dtype=object)
    new_fam = np.asarray([fams[i] for i in order], dtype=object)
    new_sex = np.asarray([sexes[i] for i in order], dtype=np.int32)
    new_father = np.full(len(order), -1, dtype=np.int32)
    new_mother = np.full(len(order), -1, dtype=np.int32)
    position = {int(old): rank for rank, old in enumerate(order.tolist())}
    for rank, old in enumerate(order.tolist()):
        fid = father_ids[old]
        mid = mother_ids[old]
        if fid is not None:
            new_father[rank] = position[id_of[fid]]
        if mid is not None:
            new_mother[rank] = position[id_of[mid]]
    # Family labels stay as strings here; Cohort construction remaps them to ints.
    return _pedigree_from_string_families(new_ids, new_father, new_mother, new_sex, new_fam)


def _pedigree_from_string_families(
    ids: np.ndarray,
    father: np.ndarray,
    mother: np.ndarray,
    sex: np.ndarray,
    families: np.ndarray,
) -> Pedigree:
    labels = {str(f): i for i, f in enumerate(dict.fromkeys(str(x) for x in families.tolist()))}
    family_id = np.asarray([labels[str(f)] for f in families.tolist()], dtype=np.int32)
    return Pedigree(
        ids=np.asarray(ids),
        father=father.astype(np.int32, copy=False),
        mother=mother.astype(np.int32, copy=False),
        sex=sex.astype(np.int32, copy=False),
        family_id=family_id,
    )


def simulate_pedigree(
    n_families: int,
    generations: int,
    mean_children: float,
    cross_rate: float,
    rng: np.random.Generator,
    min_people: int = 24,
) -> Pedigree:
    """Build a cohort of nuclear families, then ``generations`` of offspring.

    With probability ``cross_rate`` a pairing uses spouses from different
    birth families. Parents always receive lower ids than their children.
    The draw is repeated with the same generator until ``min_people`` is met
    or eight attempts have been made.
    """

    if n_families < 1:
        raise ValueError("n_families must be positive")
    best: Pedigree | None = None
    for _ in range(8):
        ped = _draw_pedigree(n_families, generations, mean_children, cross_rate, rng)
        best = ped
        if ped.n >= min_people:
            return ped
    assert best is not None
    return best


def _draw_pedigree(
    n_families: int,
    generations: int,
    mean_children: float,
    cross_rate: float,
    rng: np.random.Generator,
) -> Pedigree:
    ids: list[str] = []
    father: list[int] = []
    mother: list[int] = []
    sex: list[int] = []
    family: list[int] = []

    def add(p: int, q: int, sx: int, fam: int) -> int:
        i = len(ids)
        ids.append(f"I{i}")
        father.append(p)
        mother.append(q)
        sex.append(sx)
        family.append(fam)
        return i

    couples: list[tuple[int, int, int]] = []
    for fam in range(n_families):
        male = add(-1, -1, 1, fam)
        female = add(-1, -1, 2, fam)
        couples.append((male, female, fam))

    for _gen in range(generations):
        kids: list[tuple[int, int, int]] = []
        for male, female, fam in couples:
            n_kids = int(rng.poisson(mean_children))
            if _gen == 0:
                n_kids = max(n_kids, 1)
            for _ in range(n_kids):
                sx = 1 if rng.random() < 0.5 else 2
                # Children keep the father's birth-family label. Cross-family
                # pairing still links the two sides through the child.
                kids.append((add(male, female, sx, fam), sx, fam))
        males = [k for k in kids if k[1] == 1]
        females = [k for k in kids if k[1] == 2]
        rng.shuffle(males)
        rng.shuffle(females)
        couples = []
        for male, female in zip(males, females):
            fam = int(female[2] if rng.random() < cross_rate else male[2])
            couples.append((male[0], female[0], fam))
        if not couples:
            break

    return Pedigree(
        ids=np.asarray(ids, dtype=object),
        father=np.asarray(father, dtype=np.int32),
        mother=np.asarray(mother, dtype=np.int32),
        sex=np.asarray(sex, dtype=np.int32),
        family_id=np.asarray(family, dtype=np.int32),
    )
