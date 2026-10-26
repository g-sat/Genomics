# Pedigree LSM

`src/pedpack/pedigree_lsm.py`. Block 200 bytes, sample 40 bytes. Forty founder pairs, then one child each. The tail seals into level-0 runs. A family that sits in more than one block is copied into its own level-1 run. Old runs stay.

| Stage | Bytes written | Compaction payload | Largest family copy | Full rebuild | Runs |
| --- | ---: | ---: | ---: | ---: | --- |
| After the flush | 4800 | 0 | 0 | 4800 | 24 level-0 |
| After compacting split families | 8000 | 1920 | 120 | 4800 | 24 level-0, 16 level-1 |

Sixteen of the forty trios crossed a block boundary. Each of those copies is 120 bytes. One full rewrite of the 120 people is 4800 bytes. The log grows to 8000 bytes because the level-0 blocks are kept after the aliases move.

A child appended after its family has been compacted does not move the sealed parents. That check is in `tests/test_pedigree_lsm.py`.
