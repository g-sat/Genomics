"""Block-level measurement of a prepared layout.

One task is scheduled per touched file and placed on the live node that holds
the most replicas of the blocks that file must read. A block whose replicas
miss that node is a remote read. A block whose replicas all sit on failed
nodes is unavailable and is charged a second remote read.

This is a placement model, not a running NameNode. It exists so the analytical
formula and the learned ranker can be scored against the same rules a later
HDFS run would be judged by: bytes, remote bytes, shuffle, file count.
"""

from __future__ import annotations

from pedpack.cost import cover_group_indices
from pedpack.types import Cluster, Measurement, PreparedLayout, Query


def replica_nodes(file_id: int, block_id: int, cluster: Cluster) -> list[int]:
    """Deterministic replica set. Extra replicas occupy successive node ids."""

    n_nodes = cluster.n_nodes
    replicas = min(cluster.replication, n_nodes)
    origin = (file_id * 10007 + block_id * 97) % n_nodes
    return [(origin + step) % n_nodes for step in range(replicas)]


def measure(layout: PreparedLayout, query: Query, cluster: Cluster) -> Measurement:
    chosen, uncovered = cover_group_indices(layout, query.samples)
    if uncovered or (len(query.samples) > 0 and not chosen):
        penalty = 1.0e6 + 10.0 * max(uncovered, 1)
        return Measurement(
            bytes_read=penalty,
            remote_bytes=penalty,
            shuffle_bytes=penalty,
            seconds=penalty,
            files_touched=0,
            blocks_touched=0,
            uncovered=uncovered,
            unavailable_blocks=0,
        )

    chosen_set = set(chosen)
    q0, q1 = query.v0, query.v1
    blocks_by_file: dict[int, set[int]] = {}
    block_size = cluster.block_size
    for i, group_index in enumerate(layout.tile_group.tolist()):
        if group_index not in chosen_set:
            continue
        v0 = int(layout.tile_v0[i])
        v1 = int(layout.tile_v1[i])
        overlap_start = max(v0, q0)
        overlap_end = min(v1, q1)
        if overlap_end <= overlap_start:
            continue
        span = v1 - v0
        nbytes = float(layout.tile_nbytes[i])
        start = float(layout.tile_offset[i]) + nbytes * ((overlap_start - v0) / span)
        end = float(layout.tile_offset[i]) + nbytes * ((overlap_end - v0) / span)
        file_id = int(layout.tile_file[i])
        file_size = float(layout.file_nbytes[file_id])
        start = min(max(start, 0.0), file_size)
        end = min(max(end, start), file_size)
        if end <= start:
            end = min(file_size, start + 1.0)
        first = int(start // block_size)
        last = int((end - 1e-9) // block_size)
        touched = blocks_by_file.setdefault(file_id, set())
        touched.update(range(first, last + 1))

    failed = set(cluster.failed_nodes)
    bytes_read = 0.0
    remote_bytes = 0.0
    unavailable = 0
    blocks_touched = 0
    for file_id, blocks in blocks_by_file.items():
        file_size = float(layout.file_nbytes[file_id])
        votes: dict[int, int] = {}
        replica_map: dict[int, list[int]] = {}
        for block_id in blocks:
            replicas = replica_nodes(file_id, block_id, cluster)
            replica_map[block_id] = replicas
            for node in replicas:
                if node not in failed:
                    votes[node] = votes.get(node, 0) + 1
        task = max(votes, key=lambda node: (votes[node], -node)) if votes else None
        for block_id in blocks:
            block_start = block_id * block_size
            if block_start >= file_size:
                continue
            length = min(block_size, file_size - block_start)
            blocks_touched += 1
            bytes_read += length
            live = [node for node in replica_map[block_id] if node not in failed]
            if task is None or task not in live:
                remote_bytes += length
                if not live:
                    unavailable += 1

    chosen_sets = [layout.group_sets[i] for i in chosen]
    split_trios = 0
    for father, mother, child in query.trios:
        together = False
        for members in chosen_sets:
            if father in members and mother in members and child in members:
                together = True
                break
        if not together:
            split_trios += 1
    shuffle_bytes = split_trios * 3 * query.n_variants_span * cluster.bytes_per_genotype
    local_bytes = max(0.0, bytes_read - remote_bytes)
    seconds = (
        local_bytes / cluster.bw_local
        + remote_bytes / cluster.bw_remote
        + shuffle_bytes / cluster.bw_shuffle
        + len(blocks_by_file) * cluster.t_open
        + (bytes_read / cluster.bytes_per_genotype) / cluster.genotypes_per_second
    )
    if unavailable:
        seconds += unavailable * (block_size / cluster.bw_remote)
    return Measurement(
        bytes_read=bytes_read,
        remote_bytes=remote_bytes,
        shuffle_bytes=shuffle_bytes,
        seconds=seconds,
        files_touched=len(blocks_by_file),
        blocks_touched=blocks_touched,
        uncovered=0,
        unavailable_blocks=unavailable,
    )


def measure_many(layout: PreparedLayout, queries: list[Query], cluster: Cluster) -> list[Measurement]:
    return [measure(layout, query, cluster) for query in queries]
