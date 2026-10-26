"""The large-payload inequality matches the rounded two-layout model on sign."""

from pedpack.boundary import family_versus_variant
from pedpack.theorem import block_ratio, objective_gap
from pedpack.types import Cluster


def _cluster() -> Cluster:
    return Cluster(
        n_nodes=8,
        block_size=65536,
        replication=3,
        bytes_per_genotype=400,
        compress_base=0.45,
        kinship_alpha=0.55,
        t_open=0.001,
        storage_scan_equivalents=0.25,
    )


def test_known_regimes_keep_their_sign():
    cluster = _cluster()
    assert objective_gap(1.30, 0.05, cluster) < 0
    assert objective_gap(1.43, 0.05, cluster) > 0
    assert objective_gap(1.30, 0.60, cluster) > 0


def test_grid_disagrees_only_on_near_ties():
    cluster = _cluster()
    misses = []
    for step in range(16):
        duplication = 1.0 + 0.05 * step
        for fraction_i in range(1, 13):
            fraction = fraction_i / 20.0
            row = family_versus_variant(1200, 8000, duplication, fraction, cluster)
            limit_wins = objective_gap(duplication, fraction, cluster) < 0
            model_wins = row["family_wins"] == 1.0
            if limit_wins != model_wins:
                misses.append(row["ratio_variant_over_family"])
    assert len(misses) == 2
    assert all(abs(ratio - 1.0) < 0.01 for ratio in misses)


def test_public_window_collapses_at_a_128mb_block():
    # Phase 3 chromosome 22: 2,504 samples, 200-site window, one byte per dosage.
    # A pair is the public family. 64 KB reproduces the measured 8 blocks vs 1.
    cohort = 2504 * 200 * 1
    family = 2 * 200 * 1
    all_blocks, family_blocks, ratio = block_ratio(cohort, family, 65536)
    assert (all_blocks, family_blocks, ratio) == (8, 1, 8.0)
    big = 128 * 1024 * 1024
    all_blocks, family_blocks, ratio = block_ratio(cohort, family, big)
    assert (all_blocks, family_blocks, ratio) == (1, 1, 1.0)


def test_scaled_selective_read_is_four_blocks_at_128mb():
    # Stored byte in the scaled model is 400 * 0.45. A 30% window of 8,000
    # variants over 1,200 samples, duplication 1.30, family fraction 0.05.
    stored = 180
    window = int(0.30 * 8000)
    cohort = 1200 * window * stored
    selective = int(0.05 * 1.30 * cohort)
    all_blocks, family_blocks, ratio = block_ratio(cohort, selective, 128 * 1024 * 1024)
    assert all_blocks == 4
    assert family_blocks == 1
    assert ratio == 4.0
