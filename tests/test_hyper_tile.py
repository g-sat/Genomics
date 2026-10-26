"""Hyper-tiles on synthetic coordinates."""

from pedpack.hyper_tile import (
    choose_hyper_layout,
    cosine,
    hyper_bytes,
    phenotype_groups,
    split_bytes,
)

PEDIGREE = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
VECTORS = [
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
]
FAMILY = [0, 1, 2]
SYMPTOM = [0, 3, 6]
GENO = 100
PHENO = 20


def test_cosine_components_cut_across_families():
    assert cosine((1.0, 0.0), (0.0, 1.0)) == 0.0
    groups = phenotype_groups(VECTORS, 0.5)
    assert groups == [[0, 3, 6], [1, 4, 7], [2, 5, 8]]


def test_joint_queries_pick_the_matching_hyper_tile():
    phenotype = phenotype_groups(VECTORS, 0.5)
    pedigree_cost = (
        hyper_bytes(PEDIGREE, FAMILY, GENO, PHENO),
        hyper_bytes(PEDIGREE, SYMPTOM, GENO, PHENO),
    )
    phenotype_cost = (
        hyper_bytes(phenotype, FAMILY, GENO, PHENO),
        hyper_bytes(phenotype, SYMPTOM, GENO, PHENO),
    )
    joint_split = (
        split_bytes(PEDIGREE, phenotype, FAMILY, GENO, PHENO, "joint"),
        split_bytes(PEDIGREE, phenotype, SYMPTOM, GENO, PHENO, "joint"),
    )
    assert pedigree_cost == (360, 1080)
    assert phenotype_cost == (1080, 360)
    assert joint_split == (480, 960)
    options = {
        "pedigree": pedigree_cost,
        "phenotype": phenotype_cost,
        "split": joint_split,
    }
    assert choose_hyper_layout(0.8, options) == "pedigree"
    assert choose_hyper_layout(0.2, options) == "phenotype"
    assert choose_hyper_layout(0.5, options) == "tie"


def test_one_modality_keeps_the_files_separate():
    phenotype = phenotype_groups(VECTORS, 0.5)
    options = {
        "pedigree": (
            hyper_bytes(PEDIGREE, FAMILY, GENO, PHENO),
            hyper_bytes(PEDIGREE, SYMPTOM, GENO, PHENO),
        ),
        "phenotype": (
            hyper_bytes(phenotype, FAMILY, GENO, PHENO),
            hyper_bytes(phenotype, SYMPTOM, GENO, PHENO),
        ),
        "split": (
            split_bytes(PEDIGREE, phenotype, FAMILY, GENO, PHENO, "mendelian"),
            split_bytes(PEDIGREE, phenotype, SYMPTOM, GENO, PHENO, "differential"),
        ),
    }
    assert options["split"] == (300, 60)
    assert choose_hyper_layout(0.8, options) == "split"
    assert choose_hyper_layout(0.2, options) == "split"
    assert choose_hyper_layout(0.0, options) == "split"
