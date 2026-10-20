"""Packed Parquet layouts for the same two designs the simulator compares.

One file holds every sample. The family side packs nuclear groups into files
of ``groups_per_file`` groups, which is the pack the search actually chose,
not one object per family. A person who sits in two groups is stored twice.
Row groups are contiguous variant windows so a query can skip the rest.

PyArrow is imported inside the functions so the rest of the package does not
require it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def column_name(group: int, sample: int) -> str:
    return f"g{group}_s{sample}"


def pack_groups(groups: list[list[int]], groups_per_file: int) -> list[list[list[int]]]:
    if groups_per_file < 1:
        raise ValueError("groups_per_file must be positive")
    return [groups[start : start + groups_per_file] for start in range(0, len(groups), groups_per_file)]


def row_group_indexes(n_rows: int, row_group_size: int, v0: int, v1: int) -> list[int]:
    """Row groups whose variant span overlaps ``[v0, v1)``."""

    if row_group_size < 1:
        raise ValueError("row_group_size must be positive")
    chosen = []
    start = 0
    index = 0
    while start < n_rows:
        end = min(n_rows, start + row_group_size)
        if start < v1 and end > v0:
            chosen.append(index)
        start = end
        index += 1
    return chosen


def write_dosage_file(
    path: Path,
    columns: dict[str, np.ndarray],
    row_group_size: int,
) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    names = list(columns)
    arrays = [pa.array(np.ascontiguousarray(columns[name], dtype=np.int8), type=pa.int8()) for name in names]
    table = pa.table(dict(zip(names, arrays, strict=True)))
    pq.write_table(table, path, row_group_size=row_group_size, compression="zstd")


def read_window(
    path: Path,
    columns: list[str],
    v0: int,
    v1: int,
    row_group_size: int,
) -> tuple[int, int, int]:
    """Return compressed bytes, entry count, and row groups opened.

    The compressed-byte total is the sum of the column chunks actually named.
    The read itself is performed so the wall clock of the caller includes it.
    """

    import pyarrow.parquet as pq

    frame = pq.ParquetFile(path)
    groups = row_group_indexes(frame.metadata.num_rows, row_group_size, v0, v1)
    schema = frame.schema_arrow
    indexes = [schema.get_field_index(name) for name in columns]
    if any(index < 0 for index in indexes):
        missing = [name for name, index in zip(columns, indexes, strict=True) if index < 0]
        raise KeyError(f"missing columns in {path}: {missing[:3]}")
    compressed = 0
    for group in groups:
        row_group = frame.metadata.row_group(group)
        for index in indexes:
            compressed += row_group.column(index).total_compressed_size
    table = frame.read_row_groups(groups, columns=columns)
    return int(compressed), int(table.num_rows * table.num_columns), len(groups)
