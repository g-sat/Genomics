"""Kinship, pedigree IO, and Mendelian transmission."""

from __future__ import annotations

import numpy as np

from pedpack.genotype import (
    build_cohort,
    compound_het_pairs,
    mendelian_error_count,
)
from pedpack.grouping import agglomerate
from pedpack.pedigree import Pedigree, kinship_matrix, parse_ped, simulate_pedigree


def test_parent_child_and_sibling_kinship():
    father = np.array([-1, -1, 0, 0], dtype=np.int32)
    mother = np.array([-1, -1, 1, 1], dtype=np.int32)
    phi = kinship_matrix(father, mother)
    assert abs(phi[0, 0] - 0.5) < 1e-12
    assert abs(phi[0, 1]) < 1e-12
    assert abs(phi[2, 0] - 0.25) < 1e-12
    assert abs(phi[2, 1] - 0.25) < 1e-12
    assert abs(phi[2, 2] - 0.5) < 1e-12
    assert abs(phi[2, 3] - 0.25) < 1e-12


def test_inbred_self_kinship_exceeds_one_half():
    # 0 and 1 are unrelated founders. 2 is their child. 3 is a child of 2 and 1
    # (parent-child mating): phi(2, 1) = 0.25, so phi(3, 3) = 0.5 * (1 + 0.25).
    father = np.array([-1, -1, 0, 2], dtype=np.int32)
    mother = np.array([-1, -1, 1, 1], dtype=np.int32)
    phi = kinship_matrix(father, mother)
    assert abs(phi[3, 3] - 0.625) < 1e-12


def test_simulated_cohort_stays_mendelian_after_column_scramble():
    rng = np.random.default_rng(4)
    pedigree = simulate_pedigree(6, 2, 2.0, 0.3, rng, min_people=20)
    cohort = build_cohort(pedigree, n_variants=300, rng=rng, permute_columns=True)
    assert mendelian_error_count(cohort.genotypes, cohort.father, cohort.mother) == 0
    child = int(np.flatnonzero(cohort.father >= 0)[0])
    parent = int(cohort.father[child])
    assert abs(cohort.kinship[child, parent] - 0.25) < 1e-9
    assert abs(cohort.kinship[child, child] - 0.5) < 1e-9


def test_compound_het_counts_trans_pairs_only():
    child = np.array([1, 1, 1], dtype=np.uint8)
    mother = np.array([1, 0, 1], dtype=np.uint8)
    father = np.array([0, 1, 1], dtype=np.uint8)
    # Sites 0 and 1 form one ordered trans pair. Site 2 is heterozygous in both parents.
    assert compound_het_pairs(child, mother, father) == 1
    assert compound_het_pairs(np.array([1, 1]), np.array([1, 1]), np.array([0, 0])) == 0


def test_ped_round_trip(tmp_path):
    text = "\n".join(
        [
            "F1 A 0 0 1 1",
            "F1 B 0 0 2 1",
            "F1 C A B 1 1",
            "F2 D 0 0 1 1",
            "F2 E 0 0 2 1",
            "F2 F D E 2 1",
        ]
    )
    path = tmp_path / "tiny.ped"
    path.write_text(text, encoding="utf-8")
    pedigree = parse_ped(path)
    assert isinstance(pedigree, Pedigree)
    assert pedigree.n == 6
    phi = kinship_matrix(pedigree.father, pedigree.mother)
    # Every non-founder has kinship 1/4 with each parent.
    for child in range(pedigree.n):
        if pedigree.father[child] < 0:
            continue
        assert abs(phi[child, pedigree.father[child]] - 0.25) < 1e-12
        assert abs(phi[child, pedigree.mother[child]] - 0.25) < 1e-12


def _brute_agglomerate(kinship: np.ndarray, max_size: int, min_kinship: float) -> list[list[int]]:
    clusters = [[i] for i in range(kinship.shape[0])]
    while len(clusters) > 1:
        best = None
        best_score = min_kinship
        for a in range(len(clusters)):
            for b in range(a + 1, len(clusters)):
                if len(clusters[a]) + len(clusters[b]) > max_size:
                    continue
                block = kinship[np.ix_(clusters[a], clusters[b])]
                score = float(block.mean())
                if score > best_score:
                    best_score = score
                    best = (a, b)
        if best is None:
            break
        a, b = best
        clusters[a] = clusters[a] + clusters[b]
        del clusters[b]
    return [sorted(c) for c in clusters]


def test_agglomerate_matches_block_mean_scan():
    rng = np.random.default_rng(3)
    pedigree = simulate_pedigree(8, 2, 2.0, 0.4, rng, min_people=24)
    phi = kinship_matrix(pedigree.father, pedigree.mother)
    fast = agglomerate(phi, max_size=6, min_kinship=0.0625)
    brute = _brute_agglomerate(phi, max_size=6, min_kinship=0.0625)
    assert fast == brute
    assert sorted(i for group in fast for i in group) == list(range(pedigree.n))
    assert all(len(group) <= 6 for group in fast)
