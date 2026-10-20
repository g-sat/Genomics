"""Closed-form cost model.

The simulator is the measurement. This model is the cheap ranking function
used before a layout is measured.

The pedigree-aware formula opens the true covering groups, rounds each file's
request up to a block, and charges shuffle from the trios in the query. The
remote fraction still does not replay task placement: it rises with replication
shortfall and with the number of files opened.

The blind formula keeps group sizes and packing, and treats membership as
random. It does not apply the kinship compression discount and it does not
know which trios share a group. That is the ranker behind
``pedpack_no_pedigree_features``.
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


def trio_colocation(layout: PreparedLayout, query: Query) -> float:
    """Fraction of this query's trios that share a group.

    A trio is co-located when one sample group contains the father, the mother,
    and the child. Variant span is not part of the test: every group is tiled
    across the variant axis. With no trios in the query the rate is 1.
    """

    if not query.trios:
        return 1.0
    hits = 0
    for father, mother, child in query.trios:
        if any(father in group and mother in group and child in group for group in layout.group_sets):
            hits += 1
    return hits / len(query.trios)


def _n_samples(layout: PreparedLayout) -> int:
    highest = 0
    for group in layout.groups:
        if group.size:
            highest = max(highest, int(group.max()) + 1)
    return max(1, highest)


def _seconds_from_parts(
    bytes_read: float,
    remote_bytes: float,
    shuffle_bytes: float,
    files_touched: int,
    cluster: Cluster,
) -> float:
    local_bytes = max(0.0, bytes_read - remote_bytes)
    return (
        local_bytes / cluster.bw_local
        + remote_bytes / cluster.bw_remote
        + shuffle_bytes / cluster.bw_shuffle
        + files_touched * cluster.t_open
        + (bytes_read / cluster.bytes_per_genotype) / cluster.genotypes_per_second
    )


def _remote_fraction(files_touched: int, cluster: Cluster) -> float:
    """Replica miss rate. One file stays near the replication baseline.

    Each extra file is its own task, placed on a different preferred node, so
    the remote fraction grows with the number of files opened.
    """

    base = 0.05 * max(0.0, 1.0 - cluster.replication / max(1, cluster.n_nodes))
    spread = max(0, files_touched - 1) / max(files_touched, 1)
    return min(0.85, base + 0.25 * spread)


def _aware_bytes(layout: PreparedLayout, query: Query, cluster: Cluster) -> tuple[float, int]:
    """Payload of the true covering groups, rounded up to blocks per file."""

    chosen, uncovered = cover_group_indices(layout, query.samples)
    if uncovered or not chosen:
        return -1.0, uncovered
    chosen_set = set(chosen)
    q0, q1 = query.v0, query.v1
    per_file: dict[int, float] = {}
    for index, group_index in enumerate(layout.tile_group.tolist()):
        if group_index not in chosen_set:
            continue
        v0 = int(layout.tile_v0[index])
        v1 = int(layout.tile_v1[index])
        overlap = min(v1, q1) - max(v0, q0)
        if overlap <= 0:
            continue
        share = float(layout.tile_nbytes[index]) * (overlap / (v1 - v0))
        file_id = int(layout.tile_file[index])
        per_file[file_id] = per_file.get(file_id, 0.0) + share
    bytes_read = 0.0
    block = cluster.block_size
    for file_id, requested in per_file.items():
        file_size = float(layout.file_nbytes[file_id])
        if requested <= 0 or file_size <= 0:
            continue
        rounded = min(file_size, np.ceil(requested / block) * block)
        bytes_read += max(requested, rounded)
    return bytes_read, 0


def _blind_bytes(layout: PreparedLayout, query: Query, cluster: Cluster) -> tuple[float, float]:
    """Bytes as if group membership were a random partition of the same sizes.

    Kinship does not shrink the tile, and the identity of the queried samples
    does not change which groups are opened. Two layouts with the same group
    sizes receive the same estimate.
    """

    n_samples = _n_samples(layout)
    fraction = min(1.0, len(query.samples) / n_samples)
    bytes_read = 0.0
    groups_hit = 0.0
    for group in layout.groups:
        size = int(group.size)
        if size <= 0:
            continue
        probability = 1.0 - (1.0 - fraction) ** size
        groups_hit += probability
        bytes_read += probability * payload_bytes(
            n_samples=size,
            n_variants=query.n_variants_span,
            bytes_per_genotype=cluster.bytes_per_genotype,
            mean_kinship=0.0,
            compress_base=cluster.compress_base,
            kinship_alpha=0.0,
        )
    return bytes_read, groups_hit


def analytical_query(
    layout: PreparedLayout,
    query: Query,
    cluster: Cluster,
    pedigree_aware: bool = True,
) -> Measurement:
    """Predict one query without walking replicas.

    ``pedigree_aware`` uses the real covering groups, kinship-compressed tile
    sizes, block rounding, and this query's trio co-location. The blind model
    keeps packing and group size, and drops who is related to whom.
    """

    if pedigree_aware:
        bytes_read, uncovered = _aware_bytes(layout, query, cluster)
        if uncovered or bytes_read < 0:
            penalty = 1.0e6 + 10.0 * max(uncovered, 1)
            return Measurement(
                bytes_read=penalty,
                remote_bytes=penalty,
                shuffle_bytes=penalty,
                seconds=penalty,
                files_touched=layout.stats.n_files,
                blocks_touched=layout.stats.n_blocks,
                uncovered=max(uncovered, 0),
                unavailable_blocks=0,
            )
        coloc = trio_colocation(layout, query)
        tiles_per_file = max(1.0, layout.stats.n_tiles / max(1, layout.stats.n_files))
        chosen, _ = cover_group_indices(layout, query.samples)
        n_tiles_span = max(1, int(np.ceil(query.n_variants_span / max(1, layout.stats.width))))
        tiles_hit = max(1, len(chosen) * n_tiles_span)
        files_touched = int(min(layout.stats.n_files, max(1, round(tiles_hit / tiles_per_file))))
    else:
        bytes_read, groups_hit = _blind_bytes(layout, query, cluster)
        mean_size = max(1.0, layout.stats.mean_group_size)
        coloc = min(1.0, mean_size / _n_samples(layout))
        tiles_per_file = max(1.0, layout.stats.n_tiles / max(1, layout.stats.n_files))
        n_tiles_span = max(1, int(np.ceil(query.n_variants_span / max(1, layout.stats.width))))
        tiles_hit = max(1.0, groups_hit * n_tiles_span)
        files_touched = int(min(layout.stats.n_files, max(1, round(tiles_hit / tiles_per_file))))
        uncovered = 0

    remote_fraction = _remote_fraction(files_touched, cluster)
    remote_bytes = bytes_read * remote_fraction
    shuffle_bytes = (
        (1.0 - coloc)
        * len(query.trios)
        * query.n_variants_span
        * 3.0
        * cluster.bytes_per_genotype
    )
    return Measurement(
        bytes_read=bytes_read,
        remote_bytes=remote_bytes,
        shuffle_bytes=shuffle_bytes,
        seconds=_seconds_from_parts(bytes_read, remote_bytes, shuffle_bytes, files_touched, cluster),
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
