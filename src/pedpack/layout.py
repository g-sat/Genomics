"""Logical tiles and the files they are packed into.

A logical tile is one sample group over a contiguous variant span. Physical
files concatenate tiles. The simulator, not this module, decides which HDFS
blocks a query actually touches.
"""

from __future__ import annotations

import numpy as np

from pedpack.cost import payload_bytes
from pedpack.grouping import grouping_stats
from pedpack.types import Candidate, Cluster, Cohort, LayoutStats, PreparedLayout


def variant_windows(n_variants: int, width: int) -> list[tuple[int, int]]:
    if width < 1:
        raise ValueError("width must be positive")
    return [(start, min(start + width, n_variants)) for start in range(0, n_variants, width)]


def windows_by_bp(positions: np.ndarray, bp_width: int) -> list[tuple[int, int]]:
    """Cover ``positions`` with windows of ``bp_width`` bases.

    The last window is short when the map does not divide evenly. Every
    variant index appears in exactly one window.
    """

    if bp_width < 1:
        raise ValueError("bp_width must be positive")
    positions = np.asarray(positions)
    n = int(positions.shape[0])
    windows: list[tuple[int, int]] = []
    start = 0
    while start < n:
        end_bp = int(positions[start]) + bp_width
        end = int(np.searchsorted(positions, end_bp, side="left"))
        end = max(end, start + 1)
        end = min(end, n)
        windows.append((start, end))
        start = end
    return windows


def iter_tiles(
    n_groups: int,
    windows: list[tuple[int, int]],
    order: str,
) -> list[tuple[int, int, int]]:
    """Return ``(group_index, v0, v1)`` in physical write order."""

    if order not in {"group_major", "locus_major"}:
        raise ValueError(f"unknown pack order {order}")
    records: list[tuple[int, int, int]] = []
    if order == "group_major":
        for group_index in range(n_groups):
            for v0, v1 in windows:
                records.append((group_index, v0, v1))
    else:
        for v0, v1 in windows:
            for group_index in range(n_groups):
                records.append((group_index, v0, v1))
    return records


def _pack_files(nbytes: np.ndarray, pack: str) -> list[list[int]]:
    n = int(nbytes.shape[0])
    if n == 0:
        return []
    if pack == "naive":
        return [[i] for i in range(n)]
    if pack.startswith("count:"):
        k = int(pack.split(":", 1)[1])
        if k < 1:
            raise ValueError("count pack needs a positive tile count")
        return [list(range(i, min(i + k, n))) for i in range(0, n, k)]
    if pack.startswith("bytes:"):
        target = int(pack.split(":", 1)[1])
        if target < 1:
            raise ValueError("byte pack target must be positive")
        files: list[list[int]] = []
        cur: list[int] = []
        size = 0.0
        for i in range(n):
            length = float(nbytes[i])
            if cur and size + length > target:
                files.append(cur)
                cur = []
                size = 0.0
            cur.append(i)
            size += length
        if cur:
            files.append(cur)
        return files
    raise ValueError(f"unknown pack spec {pack}")


def group_mean_kinship(group: np.ndarray, kinship: np.ndarray) -> float:
    if group.size < 2:
        return 0.0
    block = kinship[np.ix_(group, group)]
    mask = ~np.eye(group.size, dtype=bool)
    return float(block[mask].mean())


def prepare_layout(candidate: Candidate, cohort: Cohort, cluster: Cluster) -> PreparedLayout:
    """Materialize tile byte offsets for one candidate."""

    groups = [np.asarray(g, dtype=np.int32) for g in candidate.groups]
    if not groups:
        raise ValueError("layout has no sample groups")
    windows = variant_windows(cohort.n_variants, candidate.width)
    group_phi = np.asarray(
        [group_mean_kinship(g, cohort.kinship) for g in groups], dtype=np.float64
    )
    sequence = iter_tiles(len(groups), windows, candidate.order)
    tile_group = np.asarray([g for g, _v0, _v1 in sequence], dtype=np.int32)
    tile_v0 = np.asarray([v0 for _g, v0, _v1 in sequence], dtype=np.int32)
    tile_v1 = np.asarray([v1 for _g, _v0, v1 in sequence], dtype=np.int32)
    tile_nbytes = np.zeros(len(sequence), dtype=np.float64)
    for i, group_index in enumerate(tile_group.tolist()):
        span = int(tile_v1[i] - tile_v0[i])
        tile_nbytes[i] = payload_bytes(
            n_samples=int(groups[group_index].size),
            n_variants=span,
            bytes_per_genotype=cluster.bytes_per_genotype,
            mean_kinship=float(group_phi[group_index]),
            compress_base=cluster.compress_base,
            kinship_alpha=cluster.kinship_alpha,
        )
    file_members = _pack_files(tile_nbytes, candidate.pack)
    tile_file = np.zeros(len(sequence), dtype=np.int32)
    tile_offset = np.zeros(len(sequence), dtype=np.float64)
    file_nbytes = np.zeros(len(file_members), dtype=np.float64)
    for file_id, members in enumerate(file_members):
        cursor = 0.0
        for tile_id in members:
            tile_file[tile_id] = file_id
            tile_offset[tile_id] = cursor
            cursor += float(tile_nbytes[tile_id])
        file_nbytes[file_id] = cursor

    n_blocks = 0
    for size in file_nbytes.tolist():
        if size <= 0:
            continue
        n_blocks += int(np.ceil(size / cluster.block_size))
    stats_g = grouping_stats(candidate.groups, cohort)
    primary = float(file_nbytes.sum())
    storage = primary * cluster.replication
    build = build_seconds(
        n_samples=cohort.n_samples,
        n_variants=cohort.n_variants,
        primary_bytes=primary,
        trio_colocation=stats_g["trio_colocation"],
        cluster=cluster,
    )
    sizes = np.asarray([len(g) for g in groups], dtype=np.float64)
    stats = LayoutStats(
        n_groups=len(groups),
        n_tiles=len(sequence),
        n_files=len(file_members),
        n_blocks=n_blocks,
        storage_bytes=storage,
        primary_bytes=primary,
        mean_group_size=float(sizes.mean()) if sizes.size else 0.0,
        group_size_cv=stats_g["group_size_cv"],
        trio_colocation=stats_g["trio_colocation"],
        pc_cut_rate=stats_g["pc_cut_rate"],
        mean_within_kinship=stats_g["mean_within_kinship"],
        duplication=stats_g["duplication"],
        width=candidate.width,
        order=candidate.order,
        pack=candidate.pack,
        family=candidate.family,
        build_seconds=build,
        label=candidate.label,
    )
    return PreparedLayout(
        candidate=candidate,
        groups=groups,
        group_sets=[set(g.tolist()) for g in groups],
        tile_group=tile_group,
        tile_v0=tile_v0,
        tile_v1=tile_v1,
        tile_file=tile_file,
        tile_offset=tile_offset,
        tile_nbytes=tile_nbytes,
        file_nbytes=file_nbytes,
        stats=stats,
        group_kinship=group_phi,
    )


def build_seconds(
    n_samples: int,
    n_variants: int,
    primary_bytes: float,
    trio_colocation: float,
    cluster: Cluster,
) -> float:
    """One-shot conversion cost from a sample-major source matrix.

    The source is read once. Regrouping shuffles a fraction of the raw bytes
    proportional to how many parent-child trios the destination splits. The
    primary copy is written locally and each extra replica is written at
    remote bandwidth.
    """

    raw = n_samples * n_variants * cluster.bytes_per_genotype
    regroup = raw * (1.0 - trio_colocation)
    extra = primary_bytes * max(0, cluster.replication - 1)
    return (
        raw / cluster.bw_local
        + regroup / cluster.bw_shuffle
        + primary_bytes / cluster.bw_local
        + extra / cluster.bw_remote
    )


def covers_cohort(layout: PreparedLayout, n_samples: int, n_variants: int) -> bool:
    """Every sample and every variant index is present in at least one tile."""

    seen = np.zeros(n_samples, dtype=bool)
    for group in layout.groups:
        seen[group] = True
    if not bool(seen.all()):
        return False
    if layout.tile_v0.size == 0:
        return False
    covered = np.zeros(n_variants, dtype=bool)
    for v0, v1 in zip(layout.tile_v0.tolist(), layout.tile_v1.tolist()):
        covered[v0:v1] = True
    return bool(covered.all())
