"""The two-layout tradeoff has the same sign as the measured regimes."""

from pedpack.boundary import family_versus_variant
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


def test_family_wins_when_duplication_is_moderate_and_queries_are_selective():
    row = family_versus_variant(1200, 8000, 1.30, 0.05, _cluster())
    assert row["family_wins"] == 1.0
    assert row["ratio_variant_over_family"] > 1.0


def test_family_loses_when_duplication_reaches_the_three_generation_regime():
    row = family_versus_variant(670, 8000, 1.43, 0.05, _cluster())
    assert row["family_wins"] == 0.0


def test_broad_family_queries_stop_paying_for_the_extra_copies():
    row = family_versus_variant(1200, 8000, 1.30, 0.60, _cluster())
    assert row["family_wins"] == 0.0
