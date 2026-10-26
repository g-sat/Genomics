"""Merkle checks for opaque trio blocks."""

import pytest

from pedpack.private_index import document_names_a_relationship, sample_token
from pedpack.trio_accumulator import TrioAccumulator, open_block


def _trios(n: int) -> list[tuple[int, int, int]]:
    trios = []
    child = n * 2
    for family in range(n):
        trios.append((family * 2, family * 2 + 1, child))
        child += 1
    return trios


def test_store_document_hides_tokens_and_edges():
    key = b"planner-secret"
    trios = _trios(4)
    payloads = [bytes([family]) * 120 for family in range(4)]
    tree = TrioAccumulator(key)
    tree.publish(trios, payloads)
    document = tree.storage_document()
    assert not document_names_a_relationship(document)
    for sample_id in (0, 1, 8):
        assert sample_token(key, sample_id) not in document
    assert tree.fetch(trios[1]) == payloads[1]
    assert tree.fetch(trios[1]) != payloads[0]


def test_a_tampered_path_returns_nothing():
    tree = TrioAccumulator(b"planner-secret")
    trios = _trios(8)
    payloads = [b"dosage" * 20 for _ in trios]
    tree.publish(trios, payloads)
    leaf, proof = tree.proof_for(trios[3])
    sibling, side = proof[0]
    flipped = bytes([sibling[0] ^ 0x01]) + sibling[1:]
    bad = [(flipped, side), *proof[1:]]
    with pytest.raises(PermissionError):
        open_block(tree.root, tree._blocks, leaf, bad)
    with pytest.raises(PermissionError):
        tree.fetch((0, 1, 99))


def test_a_different_key_does_not_open_the_same_root():
    trios = _trios(2)
    payloads = [b"a" * 40, b"b" * 40]
    tree = TrioAccumulator(b"planner-secret")
    tree.publish(trios, payloads)
    other = TrioAccumulator(b"other-key")
    foreign = other.leaf_for(trios[0])
    _leaf, proof = tree.proof_for(trios[0])
    with pytest.raises(PermissionError):
        open_block(tree.root, tree._blocks, foreign, proof)
