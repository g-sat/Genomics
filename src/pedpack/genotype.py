"""Mendelian genotype simulation and small inheritance checks.

Transmission is locus-independent. That is enough for a layout study: tile
geometry does not depend on linkage disequilibrium. A real codec's compression
is represented separately by a kinship term in the byte model, not by this
simulator.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from pedpack.pedigree import Pedigree, kinship_matrix, parse_ped
from pedpack.types import Cohort


def simulate_genotypes(
    pedigree: Pedigree,
    n_variants: int,
    rng: np.random.Generator,
    af_low: float = 0.01,
    af_high: float = 0.5,
    bp_spacing: int = 100,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return dosages, paternal alleles, maternal alleles, and positions.

    Dosages are in ``{0, 1, 2}``. Founders draw two Bernoulli alleles at a
    frequency sampled uniformly from ``[af_low, af_high]``. A child inherits
    one random allele from each parent.
    """

    if n_variants < 1:
        raise ValueError("n_variants must be positive")
    n = pedigree.n
    af = rng.uniform(af_low, af_high, size=n_variants)
    paternal = np.zeros((n_variants, n), dtype=np.uint8)
    maternal = np.zeros((n_variants, n), dtype=np.uint8)
    for i in range(n):
        p = int(pedigree.father[i])
        q = int(pedigree.mother[i])
        if p < 0 and q < 0:
            paternal[:, i] = rng.random(n_variants) < af
            maternal[:, i] = rng.random(n_variants) < af
            continue
        if p < 0 or q < 0:
            raise ValueError("simulation requires both parents or neither")
        take_pat = rng.random(n_variants) < 0.5
        paternal[:, i] = np.where(take_pat, paternal[:, p], maternal[:, p])
        take_mat = rng.random(n_variants) < 0.5
        maternal[:, i] = np.where(take_mat, paternal[:, q], maternal[:, q])
    dosages = paternal + maternal
    positions = (np.arange(n_variants, dtype=np.int64) * bp_spacing) + bp_spacing
    return dosages, paternal, maternal, positions


def mendelian_error_count(
    genotypes: np.ndarray,
    father: np.ndarray,
    mother: np.ndarray,
) -> int:
    """Count child/parent dosage combinations that cannot occur without error.

    The test is the usual biallelic screen: a child with dosage 2 cannot have
    a parent with dosage 0, a child with dosage 0 cannot have a parent with
    dosage 2, and a heterozygous child cannot have two homozygous parents of
    the same allele.
    """

    errors = 0
    kids = np.flatnonzero((father >= 0) & (mother >= 0))
    for child in kids.tolist():
        f = int(father[child])
        m = int(mother[child])
        gc = genotypes[:, child]
        gf = genotypes[:, f]
        gm = genotypes[:, m]
        bad = ((gc == 2) & ((gf == 0) | (gm == 0))) | ((gc == 0) & ((gf == 2) | (gm == 2)))
        bad = bad | ((gc == 1) & (gf == 0) & (gm == 0)) | ((gc == 1) & (gf == 2) & (gm == 2))
        errors += int(np.count_nonzero(bad))
    return errors


def compound_het_pairs(child: np.ndarray, mother: np.ndarray, father: np.ndarray) -> int:
    """Count ordered locus pairs in a simplified trans compound-het pattern.

    Site A is carried from the mother and site B from the father when the
    mother is heterozygous and the father is homozygous reference at A, and
    the reverse holds at B, with the child heterozygous at both. This is a
    workload predicate, not a clinical caller: there is no allele-frequency
    filter, no gene model, and no phase resolution beyond the parent dosages.
    """

    child = np.asarray(child)
    mother = np.asarray(mother)
    father = np.asarray(father)
    maternal = (mother == 1) & (father == 0) & (child == 1)
    paternal = (father == 1) & (mother == 0) & (child == 1)
    return int(np.count_nonzero(maternal) * np.count_nonzero(paternal))


def build_cohort(
    pedigree: Pedigree,
    n_variants: int,
    rng: np.random.Generator,
    bp_spacing: int = 100,
    permute_columns: bool = True,
) -> Cohort:
    """Simulate dosages and optionally scramble column order.

    Scrambling stops an index-contiguous chunking policy from accidentally
    co-locating a nuclear family just because the simulator numbered that
    family sequentially.
    """

    dosages, _paternal, _maternal, positions = simulate_genotypes(
        pedigree, n_variants, rng, bp_spacing=bp_spacing
    )
    # Kinship is computed in topological order, then the matrix is permuted
    # with the columns. The recursion cannot be rerun after the scramble: a
    # parent may sit in a higher column than their child.
    phi = kinship_matrix(pedigree.father, pedigree.mother)
    n = pedigree.n
    order = rng.permutation(n) if permute_columns else np.arange(n)
    order = np.asarray(order, dtype=np.int32)
    inverse = np.empty(n, dtype=np.int32)
    inverse[order] = np.arange(n, dtype=np.int32)

    def remap_parent(parents: np.ndarray) -> np.ndarray:
        out = np.full(n, -1, dtype=np.int32)
        for new_col, old in enumerate(order.tolist()):
            parent = int(parents[old])
            if parent >= 0:
                out[new_col] = int(inverse[parent])
        return out

    return Cohort(
        ids=np.asarray(pedigree.ids)[order],
        positions=positions,
        genotypes=np.ascontiguousarray(dosages[:, order]),
        father=remap_parent(pedigree.father),
        mother=remap_parent(pedigree.mother),
        sex=np.asarray(pedigree.sex, dtype=np.int32)[order],
        family_id=np.asarray(pedigree.family_id, dtype=np.int32)[order],
        kinship=phi[np.ix_(order, order)],
    )


def cohort_from_ped(
    path: str | Path,
    n_variants: int,
    rng: np.random.Generator,
    bp_spacing: int = 100,
) -> Cohort:
    """Parse a PED file and simulate dosages on that pedigree.

    The dosages are synthetic. The pedigree structure is the caller's. This is
    the intended route for a public pedigree extract when the genotype matrix
    itself cannot be redistributed: swap ``simulate_genotypes`` for a loader
    that returns a dosage matrix in the same individual order as ``parse_ped``
    before column permutation.
    """

    pedigree = parse_ped(path)
    return build_cohort(pedigree, n_variants, rng, bp_spacing=bp_spacing)


def load_dosage_npy(path: str | Path, pedigree: Pedigree) -> np.ndarray:
    """Load a ``(variants, samples)`` dosage matrix aligned to PED order.

    Samples must match ``pedigree.n``. Values outside ``{0, 1, 2}`` are rejected
    so a probability matrix is not silently treated as a hard call.
    """

    dosages = np.load(path)
    if dosages.ndim != 2:
        raise ValueError("dosage matrix must be two-dimensional")
    if dosages.shape[1] != pedigree.n:
        raise ValueError(
            f"dosage columns ({dosages.shape[1]}) do not match pedigree ({pedigree.n})"
        )
    if not np.isin(dosages, [0, 1, 2]).all():
        raise ValueError("dosages must be in {0, 1, 2}")
    return dosages.astype(np.uint8, copy=False)
