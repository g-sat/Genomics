"""Closed-form cost model.

The simulator is the measurement. This model is the cheap ranking function
used when a layout has not been measured yet. It is intentionally short of
the simulator in two places that show up in real block stores:

- payload bytes are proportional to the variant overlap, with no rounding up
  to the block that actually gets read;
- remote-read fraction is a function of replication and node count, not of
  which node the task was placed on.

Shuffle uses the layout-wide trio co-location rate rather than the query's
own trios. Query-specific co-location is reserved for the learned ranker.
"""

from __future__ import annotations

import numpy as np

from pedpack.types import Cluster, LayoutStats, Measurement, PreparedLayout, Query


def payload_bytes(
    n_samples: int,
    n_variants: int,
    bytes_per_genotype: float,
    mean_kinship: float,
    compress_base: float,
    kinship_alpha: float,
) -> float:
    """Modeled compressed size of one tile.

    ``compress_base`` is the codec ratio on unrelated samples. Related samples
    inside a tile reduce that ratio by ``kinship_alpha * mean_kinship``. The
    term is a stand-in for a PBWT-style codec, which does get smaller on
    related haplotypes; it is not itself a codec.
    """

    raw = float(n_samples) * float(n_variants) * float(bytes_per_genotype)
    factor = compress_base * (1.0 - kinship_alpha * mean_kinship)
    factor = min(max(factor, 0.08), 1.0)
    return raw * factor


def cover_group_indices(layout: PreparedLayout, samples: tuple[int, ...] | np.ndarray) -> tuple[list[int], int]:
    """Greedy set cover of ``samples`` by this layout's groups.

    The same planner is used by the analytical model and the simulator, so
    they disagree about bytes and placement, not about which groups are chosen.
    """

    remaining = {int(s) for s in np.asarray(samples).tolist()}
    if not remaining:
        return [], 0
    chosen: list[int] = []
    while remaining:
        best_i = -1
        best_gain = 0
        for i, members in enumerate(layout.group_sets):
            gain = len(members & remaining)
            if gain > best_gain:
                best_gain = gain
                best_i = i
        if best_gain == 0:
            break
        chosen.append(best_i)
        remaining -= layout.group_sets[best_i]
    return chosen, len(remaining)


def analytical_query(layout: PreparedLayout, query: Query, cluster: Cluster) -> Measurement:
    """Predict one query without walking block replicas."""

    chosen, uncovered = cover_group_indices(layout, query.samples)
    if uncovered or not chosen:
        penalty = 1.0e6 + 10.0 * uncovered
        return Measurement(
            bytes_read=penalty,
            remote_bytes=penalty,
            shuffle_bytes=penalty,
            seconds=penalty,
            files_touched=layout.stats.n_files,
            blocks_touched=layout.stats.n_blocks,
            uncovered=uncovered,
            unavailable_blocks=0,
        )

    q0, q1 = query.v0, query.v1
    bytes_read = 0.0
    tiles_hit = 0
    groups_hit = np.zeros(layout.stats.n_groups, dtype=bool)
    for group_index in chosen:
        groups_hit[group_index] = True
    for i, group_index in enumerate(layout.tile_group.tolist()):
        if not groups_hit[group_index]:
            continue
        v0 = int(layout.tile_v0[i])
        v1 = int(layout.tile_v1[i])
        overlap = min(v1, q1) - max(v0, q0)
        if overlap <= 0:
            continue
        span = v1 - v0
        bytes_read += float(layout.tile_nbytes[i]) * (overlap / span)
        tiles_hit += 1

    remote_fraction = 0.05 * max(0.0, 1.0 - cluster.replication / cluster.n_nodes)
    remote_bytes = bytes_read * remote_fraction
    # Layout-wide co-location, not this query's trios. That gap is deliberate:
    # the learned ranker is the component that sees query-specific co-location.
    shuffle_bytes = (
        (1.0 - layout.stats.trio_colocation)
        * len(query.trios)
        * query.n_variants_span
        * 3.0
        * cluster.bytes_per_genotype
    )
    tiles_per_file = max(1.0, layout.stats.n_tiles / max(1, layout.stats.n_files))
    files_touched = int(min(layout.stats.n_files, max(1, round(tiles_hit / tiles_per_file))))
    local_bytes = bytes_read - remote_bytes
    seconds = (
        local_bytes / cluster.bw_local
        + remote_bytes / cluster.bw_remote
        + shuffle_bytes / cluster.bw_shuffle
        + files_touched * cluster.t_open
        + (bytes_read / cluster.bytes_per_genotype) / cluster.genotypes_per_second
    )
    return Measurement(
        bytes_read=bytes_read,
        remote_bytes=remote_bytes,
        shuffle_bytes=shuffle_bytes,
        seconds=seconds,
        files_touched=files_touched,
        blocks_touched=0,
        uncovered=0,
        unavailable_blocks=0,
    )


def storage_weight(reference_scan_seconds: float, reference_storage_bytes: float, equivalents: float) -> float:
    """Seconds of objective charged per replicated byte."""

    if reference_storage_bytes <= 0:
        return 0.0
    return equivalents * reference_scan_seconds / reference_storage_bytes


def objective(stats: LayoutStats, query_seconds: float, weight: float) -> float:
    """Build plus priced storage plus execution, all in seconds."""

    return stats.build_seconds + stats.storage_bytes * weight + query_seconds


def steady_state(stats: LayoutStats, mean_query_seconds: float, weight: float) -> float:
    """Objective after build has been paid, per query."""

    return stats.storage_bytes * weight + mean_query_seconds
