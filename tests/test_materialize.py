"""Packed files match the modeled tile sizes, and a family query skips other groups."""

from __future__ import annotations

import numpy as np

from pedpack.genotype import build_cohort
from pedpack.grouping import nuclear_groups
from pedpack.layout import prepare_layout
from pedpack.materialize import read_plan, write_layout
from pedpack.pedigree import simulate_pedigree
from pedpack.types import Candidate, Cluster, Query


def test_family_query_reads_less_than_the_whole_layout(tmp_path):
    rng = np.random.default_rng(1)
    pedigree = simulate_pedigree(4, 1, 2.0, 0.0, rng, min_people=12)
    cohort = build_cohort(pedigree, 40, rng)
    cluster = Cluster(bytes_per_genotype=32, block_size=1024, compress_base=1.0, kinship_alpha=0.0)
    nuclear = nuclear_groups(cohort)
    layout = prepare_layout(
        Candidate("f", "family", tuple(tuple(g) for g in nuclear), 40, "naive", "group_major"),
        cohort,
        cluster,
    )
    manifest = write_layout(layout, cohort, tmp_path)
    written = sum((tmp_path / record["path"]).stat().st_size for record in manifest["files"])
    assert written == manifest["primary_bytes"]
    family = nuclear[0]
    query = Query("one", "family_selective", tuple(family), 0, 40, ())
    plan = read_plan(manifest, query)
    assert sum(item["length"] for item in plan) < manifest["primary_bytes"]
