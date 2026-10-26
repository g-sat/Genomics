"""Pointer stitch against rounded block reads."""

from pedpack.pointer_resolve import PointerPool


def test_a_shared_parent_is_stored_once():
    pool = PointerPool(block_bytes=256, sample_bytes=100)
    pool.add_sample(0)
    pool.add_sample(1)
    pool.add_sample(2, father=0, mother=1)
    pool.add_sample(3, father=0, mother=1)
    assert pool.stored_bytes() == 400
    assert pool.parents[2] == (0, 1)
    assert pool.stitch_bytes([0, 1, 2, 0, 1, 3]) == 400
    assert pool.materialized_bytes([(0, 1, 2), (0, 1, 3)]) == 600


def test_stitch_skips_the_rest_of_the_block():
    pool = PointerPool(block_bytes=256, sample_bytes=100)
    for sample_id in range(4):
        pool.add_sample(sample_id)
    assert pool.rounded_read([0, 1, 2]) == 512
    assert pool.stitch_bytes([0, 1, 2]) == 300


def test_scattered_columns_open_one_block_each():
    pool = PointerPool(block_bytes=200, sample_bytes=100)
    pool.add_sample(0)
    pool.add_sample(10)
    pool.add_sample(1)
    pool.add_sample(11)
    pool.add_sample(2, father=0, mother=1)
    assert pool.rounded_read([0, 1, 2]) == 600
    assert pool.stitch_bytes([0, 1, 2]) == 300
