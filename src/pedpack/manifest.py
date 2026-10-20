"""JSON description of a packed layout.

This is the handoff a storage writer would consume: one record per physical
file, with tile offsets in bytes. It does not open a NameNode.
"""

from __future__ import annotations

import numpy as np

from pedpack.types import Cluster, PreparedLayout


def layout_manifest(layout: PreparedLayout, cluster: Cluster) -> dict:
    files = []
    for file_id, size in enumerate(layout.file_nbytes.tolist()):
        indexes = np.nonzero(layout.tile_file == file_id)[0]
        tiles = []
        for index in indexes.tolist():
            group = layout.groups[int(layout.tile_group[index])]
            tiles.append(
                {
                    "group_index": int(layout.tile_group[index]),
                    "samples": [int(s) for s in group.tolist()],
                    "variant_start": int(layout.tile_v0[index]),
                    "variant_end": int(layout.tile_v1[index]),
                    "offset": float(layout.tile_offset[index]),
                    "nbytes": float(layout.tile_nbytes[index]),
                }
            )
        files.append({"file_id": int(file_id), "nbytes": float(size), "tiles": tiles})
    return {
        "label": layout.stats.label,
        "family": layout.stats.family,
        "block_size": int(cluster.block_size),
        "replication": int(cluster.replication),
        "n_files": int(layout.stats.n_files),
        "n_blocks": int(layout.stats.n_blocks),
        "files": files,
    }
