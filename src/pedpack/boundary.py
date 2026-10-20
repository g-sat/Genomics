"""When a nuclear-family layout beats an all-samples layout.

The full ranker still scores every candidate. This function is the
two-layout slice of that trade-off: duplication against the fraction of
families a selective query names. It uses the same block rounding, remote
fraction, build model, and storage price as the rest of the cost model.
It does not apply the kinship compression discount, so a family win here
is not a gift from the codec term.

A query that names a fraction ``f`` of the nuclear groups reads about
``f`` of the family layout and all of the all-samples layout. Extended
queries are treated as naming twice as many families, capped at the whole
cohort. Both layouts co-locate trios, so neither pays shuffle.
"""

from __future__ import annotations

import math

from pedpack.cost import _seconds_from_parts, storage_weight
from pedpack.layout import build_seconds
from pedpack.types import Cluster


def _rounded(total: float, n_files: int, block: int) -> float:
    files = max(1, int(n_files))
    each = total / files
    per = math.ceil(each / block) * block if each > 0 else 0.0
    return per * files


def _query_seconds(nbytes: float, files_touched: int, cluster: Cluster, remote_fraction: float) -> float:
    remote = remote_fraction * nbytes
    return _seconds_from_parts(nbytes, remote, 0.0, files_touched, cluster)


def family_versus_variant(
    n_samples: int,
    n_variants: int,
    duplication: float,
    family_fraction: float,
    cluster: Cluster,
    window_fraction: float = 0.30,
    mean_group_size: float = 4.0,
    groups_per_file: int = 32,
    horizon: float = 32.0,
    scan_share: float = 0.35,
    family_share: float = 0.40,
    extended_share: float = 0.25,
    storage_scan_equivalents: float = 0.25,
) -> dict[str, float]:
    """Return both objectives and the variant/family ratio.

    A ratio above 1 means the all-samples layout costs more, so the family
    layout wins.
    """

    if n_samples < 1 or n_variants < 1:
        raise ValueError("cohort dimensions must be positive")
    if duplication < 1:
        raise ValueError("duplication is stored copies per sample and is at least 1")
    fraction = min(1.0, max(0.0, float(family_fraction)))
    window = min(1.0, max(1.0 / n_variants, float(window_fraction)))
    payload = cluster.bytes_per_genotype * cluster.compress_base
    raw = n_samples * n_variants * payload
    n_groups = max(1, math.ceil(duplication * n_samples / mean_group_size))
    family_files = max(1, math.ceil(n_groups / groups_per_file))
    variant_primary = _rounded(raw, 1, cluster.block_size)
    family_primary = _rounded(raw * duplication, family_files, cluster.block_size)
    variant_storage = variant_primary * cluster.replication
    family_storage = family_primary * cluster.replication

    def opened_files(n_groups_touched: int) -> int:
        """How many packed files a scattered set of groups lands in."""

        if n_groups_touched >= n_groups:
            return family_files
        per = groups_per_file
        # 1 - ((G - per) / G) ^ k, times the file count.
        stay = max(0.0, (n_groups - per) / n_groups)
        miss = stay ** n_groups_touched
        return max(1, min(family_files, math.ceil(family_files * (1.0 - miss))))

    def group_bytes(share: float) -> tuple[float, int, int]:
        touched_groups = max(1, math.ceil(share * n_groups))
        each = mean_group_size * (window * n_variants) * payload
        each = math.ceil(each / cluster.block_size) * cluster.block_size
        return touched_groups * each, opened_files(touched_groups), touched_groups

    # On the scaled simulator both the one-file layout and the family layout
    # read with a remote fraction near 0.6. The analytical floor would hand
    # the one-file layout a placement win the measurement does not.
    remote_fraction = 0.60
    variant_scan = _query_seconds(variant_primary, 1, cluster, remote_fraction)
    family_scan = _query_seconds(family_primary, family_files, cluster, remote_fraction)
    variant_selective = _query_seconds(variant_primary * window, 1, cluster, remote_fraction)
    family_selective_bytes, family_selective_files, _groups = group_bytes(fraction)
    extended = min(1.0, fraction * 2.0)
    family_extended_bytes, family_extended_files, _ext_groups = group_bytes(extended)

    variant_query = (
        scan_share * variant_scan
        + family_share * variant_selective
        + extended_share * variant_selective
    )
    family_query = (
        scan_share * family_scan
        + family_share * _query_seconds(family_selective_bytes, family_selective_files, cluster, remote_fraction)
        + extended_share * _query_seconds(family_extended_bytes, family_extended_files, cluster, remote_fraction)
    )
    variant_build = build_seconds(n_samples, n_variants, variant_primary, 1.0, cluster)
    family_build = build_seconds(n_samples, n_variants, family_primary, 1.0, cluster)
    weight = storage_weight(variant_scan, variant_storage, storage_scan_equivalents)
    variant_objective = variant_build + weight * variant_storage + horizon * variant_query
    family_objective = family_build + weight * family_storage + horizon * family_query
    return {
        "variant_objective": variant_objective,
        "family_objective": family_objective,
        "ratio_variant_over_family": variant_objective / family_objective,
        "family_wins": 1.0 if family_objective < variant_objective else 0.0,
        "duplication": float(duplication),
        "family_fraction": fraction,
        "variant_build": variant_build,
        "family_build": family_build,
        "variant_storage_seconds": weight * variant_storage,
        "family_storage_seconds": weight * family_storage,
        "variant_query_seconds": horizon * variant_query,
        "family_query_seconds": horizon * family_query,
        "variant_scan": variant_scan,
        "family_scan": family_scan,
    }
