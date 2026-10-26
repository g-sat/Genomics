"""Family compaction writes one family, not the cohort."""

from pedpack.pedigree_lsm import PedigreeLSM


def _stream(n_families: int) -> PedigreeLSM:
    tree = PedigreeLSM(block_bytes=200, sample_bytes=40)
    child = n_families * 2
    for family in range(n_families):
        father = family * 2
        mother = father + 1
        tree.add_sample(father)
        tree.add_sample(mother)
        tree.add_sample(child, father=father, mother=mother)
        child += 1
    return tree


def test_a_new_child_does_not_copy_parents():
    tree = PedigreeLSM(block_bytes=200, sample_bytes=40)
    tree.add_sample(0)
    tree.add_sample(1)
    tree.add_sample(2, father=0, mother=1)
    assert tree.compaction_payload == 0
    assert tree.family_blocks(2) == {(-1, 0)}


def test_compaction_moves_only_the_family():
    tree = PedigreeLSM(block_bytes=200, sample_bytes=40)
    tree.add_sample(0)
    tree.add_sample(1)
    tree.add_sample(10)
    tree.add_sample(11)
    tree.add_sample(12)
    tree.add_sample(2, father=0, mother=1)
    tree.flush()
    assert len(tree.family_blocks(2)) == 2
    stranger = tree.alias[10]
    tree.compact_family(2)
    assert tree.compaction_payload == 3 * 40
    assert tree.largest_compaction == 3 * 40
    assert tree.largest_compaction < tree.full_rebuild_bytes()
    assert tree.alias[10] == stranger
    assert tree.family_blocks(2) == {tree.alias[0]}
    assert tree.runs[tree.alias[0][0]].level == 1
    assert tree.level_counts()[0] >= 1
    assert tree.level_counts()[1] == 1


def test_a_later_child_appends_beside_a_sealed_family():
    tree = _stream(1)
    tree.flush()
    tree.compact_family(2)
    sealed = tree.alias[0]
    tree.add_sample(3, father=0, mother=1)
    assert tree.alias[0] == sealed
    assert tree.compaction_payload == 3 * 40
    assert 3 in tree._tail
