"""Packed Parquet reads the groups a query names."""

from __future__ import annotations

import numpy as np

from pedpack.columnar import column_name, pack_groups, read_window, row_group_indexes, write_dosage_file


def test_row_groups_cover_the_window_only():
    assert row_group_indexes(100, 40, 0, 40) == [0]
    assert row_group_indexes(100, 40, 39, 41) == [0, 1]
    assert row_group_indexes(100, 40, 80, 100) == [2]


def test_family_window_is_smaller_than_the_full_matrix(tmp_path):
    rng = np.random.default_rng(0)
    genotypes = rng.integers(0, 3, size=(80, 12), dtype=np.int8)
    groups = [[0, 1, 2, 3], [4, 5, 6], [2, 7, 8], [9, 10, 11]]
    variant = tmp_path / "variant.parquet"
    write_dosage_file(variant, {str(i): genotypes[:, i] for i in range(12)}, row_group_size=20)
    packed = pack_groups(groups, groups_per_file=2)
    family = tmp_path / "part-0.parquet"
    columns = {
        column_name(index, sample): genotypes[:, sample]
        for index, group in enumerate(packed[0])
        for sample in group
    }
    write_dosage_file(family, columns, row_group_size=20)
    full_bytes, _, full_groups = read_window(variant, [str(i) for i in range(12)], 0, 80, 20)
    part_bytes, _, part_groups = read_window(
        family,
        [column_name(0, sample) for sample in groups[0]],
        20,
        40,
        20,
    )
    assert full_groups == 4
    assert part_groups == 1
    assert part_bytes < full_bytes
