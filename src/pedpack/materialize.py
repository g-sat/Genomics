"""Write a prepared layout as real files a distributed store can hold.

Each tile is padded to the byte size the cost model assigned. The leading
bytes are the dosage values for that tile, one byte per genotype, in variant-
major order. The rest of the tile is zero. There is no codec: the padding is
the modeled record width, so a range read asks for the same number of bytes
the simulator charged.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pedpack.cost import cover_group_indices
from pedpack.types import Cohort, PreparedLayout, Query


def integer_tile_sizes(layout: PreparedLayout) -> np.ndarray:
    sizes = np.rint(layout.tile_nbytes).astype(np.int64)
    return np.maximum(sizes, 1)


def _tile_payload(cohort: Cohort, samples: np.ndarray, v0: int, v1: int, nbytes: int) -> bytes:
    dosages = cohort.genotypes[v0:v1][:, samples].reshape(-1)
    raw = np.ascontiguousarray(dosages, dtype=np.uint8).tobytes()
    if len(raw) >= nbytes:
        return raw[:nbytes]
    return raw + bytes(nbytes - len(raw))


def write_layout(layout: PreparedLayout, cohort: Cohort, directory: str | Path) -> dict:
    """Write one file per packed HDFS file. Returns integer offsets."""

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    sizes = integer_tile_sizes(layout)
    n_files = int(layout.stats.n_files)
    file_records = []
    for file_id in range(n_files):
        indexes = np.nonzero(layout.tile_file == file_id)[0]
        indexes = indexes[np.argsort(layout.tile_offset[indexes], kind="mergesort")]
        chunks: list[bytes] = []
        tiles = []
        cursor = 0
        for index in indexes.tolist():
            group = layout.groups[int(layout.tile_group[index])]
            v0 = int(layout.tile_v0[index])
            v1 = int(layout.tile_v1[index])
            nbytes = int(sizes[index])
            chunks.append(_tile_payload(cohort, group, v0, v1, nbytes))
            tiles.append(
                {
                    "tile": int(index),
                    "group_index": int(layout.tile_group[index]),
                    "samples": [int(s) for s in group.tolist()],
                    "variant_start": v0,
                    "variant_end": v1,
                    "offset": cursor,
                    "nbytes": nbytes,
                }
            )
            cursor += nbytes
        blob = b"".join(chunks)
        path = directory / f"part-{file_id:05d}"
        path.write_bytes(blob)
        file_records.append({"file_id": file_id, "path": path.name, "nbytes": cursor, "tiles": tiles})
    manifest = {
        "label": layout.stats.label,
        "family": layout.stats.family,
        "width": int(layout.stats.width),
        "pack": layout.stats.pack,
        "order": layout.stats.order,
        "n_files": n_files,
        "n_blocks_model": int(layout.stats.n_blocks),
        "storage_bytes_model": float(layout.stats.storage_bytes),
        "primary_bytes": int(sum(record["nbytes"] for record in file_records)),
        "duplication": float(layout.stats.duplication),
        "trio_colocation": float(layout.stats.trio_colocation),
        "files": file_records,
    }
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return manifest


def read_plan(manifest: dict, query: Query) -> list[dict]:
    """Byte ranges a query should open. Ranges follow the variant overlap inside each tile."""

    class _Groups:
        group_sets = [set(tile["samples"]) for record in manifest["files"] for tile in record["tiles"]]

    # Cover is over groups, which repeat on every tile. Use the first tile of each group.
    group_samples: dict[int, set[int]] = {}
    for record in manifest["files"]:
        for tile in record["tiles"]:
            group_samples.setdefault(tile["group_index"], set(tile["samples"]))
    sets = [group_samples[i] for i in range(len(group_samples))]

    class _Layout:
        group_sets = sets

    chosen, uncovered = cover_group_indices(_Layout(), query.samples)
    if uncovered or not chosen:
        raise ValueError(f"query {query.name} is not covered ({uncovered} samples)")
    chosen_set = set(chosen)
    ranges = []
    for record in manifest["files"]:
        for tile in record["tiles"]:
            if tile["group_index"] not in chosen_set:
                continue
            v0, v1 = tile["variant_start"], tile["variant_end"]
            overlap = min(v1, query.v1) - max(v0, query.v0)
            if overlap <= 0:
                continue
            span = v1 - v0
            length = max(1, int(round(tile["nbytes"] * (overlap / span))))
            start_frac = (max(v0, query.v0) - v0) / span
            offset = tile["offset"] + int(tile["nbytes"] * start_frac)
            length = min(length, tile["nbytes"] - int(tile["nbytes"] * start_frac))
            ranges.append(
                {
                    "query": query.name,
                    "file": record["path"],
                    "offset": int(offset),
                    "length": int(max(1, length)),
                }
            )
    return ranges
