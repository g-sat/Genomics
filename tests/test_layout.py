"""Layout geometry, packing, and the gap between the formula and the simulator."""

from __future__ import annotations

import numpy as np

from pedpack.cost import analytical_query
from pedpack.genotype import build_cohort
from pedpack.grouping import generic_chunks, grouping_stats, nuclear_groups
from pedpack.layout import covers_cohort, prepare_layout, windows_by_bp
from pedpack.manifest import layout_manifest
from pedpack.pedigree import simulate_pedigree, trios
from pedpack.simulator import measure, replica_nodes
from pedpack.types import Candidate, Cluster, Cohort, Query


def _founders(n: int, m: int) -> Cohort:
    return Cohort(
        ids=np.arange(n),
        positions=np.arange(1, m + 1, dtype=np.int64) * 100,
        genotypes=np.zeros((m, n), dtype=np.uint8),
        father=np.full(n, -1, dtype=np.int32),
        mother=np.full(n, -1, dtype=np.int32),
        sex=np.ones(n, dtype=np.int32),
        family_id=np.zeros(n, dtype=np.int32),
        kinship=np.eye(n, dtype=np.float64) * 0.5,
    )


def test_bp_windows_cover_every_variant_once():
    positions = np.array([100, 150, 400, 1000, 1001], dtype=np.int64)
    windows = windows_by_bp(positions, bp_width=300)
    covered = []
    for start, end in windows:
        covered.extend(range(start, end))
    assert covered == [0, 1, 2, 3, 4]


def test_packed_narrow_query_reads_the_shared_block():
    cohort = _founders(10, 100)
    cluster = Cluster(
        n_nodes=4,
        block_size=1000,
        replication=1,
        bytes_per_genotype=1.0,
        compress_base=1.0,
        kinship_alpha=0.0,
        t_open=0.0,
    )
    candidate = Candidate(
        label="pack",
        family="variant",
        groups=(tuple(range(10)),),
        width=50,
        pack="count:2",
        order="group_major",
    )
    layout = prepare_layout(candidate, cohort, cluster)
    assert layout.stats.n_files == 1
    assert covers_cohort(layout, 10, 100)
    query = Query("narrow", "fit", tuple(range(10)), 0, 50, ())
    simulated = measure(layout, query, cluster)
    formula = analytical_query(layout, query, cluster)
    # The queried half of the file shares its HDFS block with the other tile.
    assert simulated.bytes_read == 1000
    assert formula.bytes_read == 1000


def test_naive_pack_creates_one_file_per_tile():
    cohort = _founders(4, 40)
    cluster = Cluster(bytes_per_genotype=1.0, compress_base=1.0, kinship_alpha=0.0, block_size=64)
    tiles = Candidate(
        label="naive",
        family="variant",
        groups=(tuple(range(4)),),
        width=10,
        pack="naive",
        order="group_major",
    )
    packed = Candidate(
        label="packed",
        family="variant",
        groups=(tuple(range(4)),),
        width=10,
        pack="count:4",
        order="group_major",
    )
    naive = prepare_layout(tiles, cohort, cluster)
    dense = prepare_layout(packed, cohort, cluster)
    assert naive.stats.n_files == naive.stats.n_tiles
    assert dense.stats.n_files < naive.stats.n_files


def test_nuclear_groups_cut_a_selective_read():
    rng = np.random.default_rng(1)
    pedigree = simulate_pedigree(6, 1, 2.2, 0.2, rng, min_people=16)
    cohort = build_cohort(pedigree, 200, rng)
    cluster = Cluster(
        block_size=4096,
        bytes_per_genotype=8.0,
        compress_base=1.0,
        kinship_alpha=0.0,
        replication=1,
        n_nodes=4,
    )
    nuclear = nuclear_groups(cohort)
    stats = grouping_stats(nuclear, cohort)
    assert stats["trio_colocation"] == 1.0
    generic = grouping_stats(generic_chunks(cohort.n_samples, 4), cohort)
    assert generic["trio_colocation"] < 0.5

    family = nuclear[0]
    query = Query("one-family", "selective", tuple(family), 0, cohort.n_variants, ())
    variant = prepare_layout(
        Candidate("v", "variant", (tuple(range(cohort.n_samples)),), cohort.n_variants, "naive", "group_major"),
        cohort,
        cluster,
    )
    family_layout = prepare_layout(
        Candidate("f", "family", tuple(tuple(g) for g in nuclear), cohort.n_variants, "naive", "group_major"),
        cohort,
        cluster,
    )
    variant_bytes = measure(variant, query, cluster).bytes_read
    family_bytes = measure(family_layout, query, cluster).bytes_read
    assert family_bytes < variant_bytes * 0.5


def test_manifest_lists_every_file():
    cohort = _founders(4, 40)
    cluster = Cluster(bytes_per_genotype=1.0, compress_base=1.0, kinship_alpha=0.0, block_size=64)
    layout = prepare_layout(
        Candidate("naive", "variant", (tuple(range(4)),), 10, "naive", "group_major"),
        cohort,
        cluster,
    )
    manifest = layout_manifest(layout, cluster)
    assert manifest["n_files"] == layout.stats.n_files
    assert len(manifest["files"]) == layout.stats.n_files
    for record, size in zip(manifest["files"], layout.file_nbytes.tolist()):
        assert abs(sum(tile["nbytes"] for tile in record["tiles"]) - size) < 1e-6


def test_failed_replica_is_a_remote_read():
    cohort = _founders(4, 20)
    cluster = Cluster(
        n_nodes=4,
        replication=1,
        block_size=1024,
        bytes_per_genotype=1.0,
        compress_base=1.0,
        kinship_alpha=0.0,
    )
    layout = prepare_layout(
        Candidate("v", "variant", (tuple(range(4)),), 20, "naive", "group_major"),
        cohort,
        cluster,
    )
    holder = replica_nodes(0, 0, cluster)[0]
    failed = cluster.with_updates(failed_nodes=(holder,))
    query = Query("all", "scan", tuple(range(4)), 0, 20, ())
    healthy = measure(layout, query, cluster)
    broken = measure(layout, query, failed)
    assert healthy.remote_fraction == 0.0
    assert broken.remote_fraction == 1.0
    assert broken.unavailable_blocks >= 1


def test_blind_formula_misses_a_colocated_trio():
    rng = np.random.default_rng(2)
    pedigree = simulate_pedigree(3, 1, 2.0, 0.0, rng, min_people=8)
    cohort = build_cohort(pedigree, 40, rng)
    cluster = Cluster(bytes_per_genotype=50.0, block_size=1024, compress_base=0.5, kinship_alpha=0.4)
    nuclear = nuclear_groups(cohort)
    layout = prepare_layout(
        Candidate("f", "family", tuple(tuple(g) for g in nuclear), 40, "naive", "group_major"),
        cohort,
        cluster,
    )
    family = set(nuclear[0])
    inside = [trio for trio in trios(cohort.father, cohort.mother) if set(trio) <= family]
    assert inside
    query = Query("q", "family_selective", tuple(sorted(family)), 0, 20, (inside[0],))
    aware = analytical_query(layout, query, cluster, pedigree_aware=True)
    blind = analytical_query(layout, query, cluster, pedigree_aware=False)
    assert aware.shuffle_bytes == 0.0
    assert blind.shuffle_bytes > 0.0
    assert aware.bytes_read < blind.bytes_read
