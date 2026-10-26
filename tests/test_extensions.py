"""Incremental aliasing, a token index, and the phenotype cut."""

from pedpack.incremental import AppendOnlyLayout
from pedpack.multimodal import choose_cut, group_bytes, pedigree_threshold
from pedpack.private_index import TokenIndex, document_names_a_relationship, sample_token


def test_new_child_does_not_rewrite_parents():
    layout = AppendOnlyLayout(block_bytes=200, sample_bytes=40)
    layout.add_sample(0)
    layout.add_sample(1)
    layout.add_sample(2, father=0, mother=1)
    assert layout.rewritten_bytes == 0
    assert layout.trio_blocks(2) == {-1}
    assert layout.alias[0] == layout.alias[1] == layout.alias[2]


def test_sealing_a_family_copies_only_that_trio():
    layout = AppendOnlyLayout(block_bytes=200, sample_bytes=40)
    for sample_id in range(6):
        if sample_id >= 2 and sample_id % 3 == 2:
            layout.add_sample(sample_id, father=sample_id - 2, mother=sample_id - 1)
        else:
            layout.add_sample(sample_id)
    layout.seal_family(2)
    assert layout.rewritten_bytes == 3 * 40
    assert layout.trio_blocks(2) == {layout.alias[0]}
    assert layout.alias[0] == layout.alias[1] == layout.alias[2]
    assert layout.full_rebuild_bytes() > layout.rewritten_bytes


def test_storage_document_has_no_pedigree_words():
    key = b"planner-secret"
    index = TokenIndex(key)
    index.publish(0, [0, 1])
    index.publish(1, [2])
    document = index.storage_document()
    assert not document_names_a_relationship(document)
    tokens = index.plan([0, 1, 2])
    assert index.blocks_for(tokens) == [0, 0, 1]
    assert sample_token(key, 0) == tokens[0]
    assert sample_token(b"other-key", 0) != tokens[0]


def test_cut_follows_the_hereditary_fraction():
    pedigree = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    phenotype = [[0, 3, 6], [1, 4, 7], [2, 5, 8]]
    sample_bytes = 100
    h_ped = group_bytes(pedigree, [0, 1, 2], sample_bytes)
    h_phe = group_bytes(phenotype, [0, 1, 2], sample_bytes)
    p_ped = group_bytes(pedigree, [0, 3, 6], sample_bytes)
    p_phe = group_bytes(phenotype, [0, 3, 6], sample_bytes)
    assert (h_ped, h_phe, p_ped, p_phe) == (300, 900, 900, 300)
    threshold = pedigree_threshold(h_ped, h_phe, p_ped, p_phe)
    assert threshold == 0.5
    assert choose_cut(0.8, h_ped, h_phe, p_ped, p_phe) == "pedigree"
    assert choose_cut(0.2, h_ped, h_phe, p_ped, p_phe) == "phenotype"
    assert choose_cut(0.5, h_ped, h_phe, p_ped, p_phe) == "tie"
