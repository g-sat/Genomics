# Large-payload condition

`src/pedpack/theorem.py` drops block rounding and file opens. Both layouts keep remote fraction 0.60. Kinship and shuffle are zero. Storage of the all-samples layout is priced at 0.25 of its scan. Horizon is 32. The family layout is cheaper when

`(D - 1) * (c_write + 0.25 * c_byte) + 32 * c_byte * Q < 0`

with `c_byte = 9.625e-9` and `c_write = 3e-8` on the paper's bandwidths.

| D | f | Gap per payload byte | Rounded variant/family ratio |
| ---: | ---: | ---: | ---: |
| 1.30 | 0.05 | −1.259e-8 | 1.061 |
| 1.43 | 0.05 | +6.175e-9 | 0.970 |
| 1.30 | 0.60 | +4.086e-8 | 0.833 |

A negative gap means the family layout wins. On D from 1.00 to 1.75 step 0.05 and f from 0.05 to 0.60 step 0.05, at 1,200 samples and 8,000 variants, the sign matches `family_versus_variant` on 190 of 192 cells. The misses are ratios 0.995 and 0.998.

Block counts use `block_ratio`. Stored byte 180 is `400 * 0.45`.

| Read | Block | Blocks, all samples | Blocks, comparison | Ratio |
| --- | ---: | ---: | ---: | ---: |
| 2,504 × 200 × 1 byte, against a pair | 64 KB | 8 | 1 | 8 |
| same window | 128 MB | 1 | 1 | 1 |
| 1,200 × 2,400 × 180 bytes, against 5% at D=1.30 | 128 MB | 4 | 1 | 4 |
| 1,200 × 100,000 × 180 bytes, same fraction | 128 MB | 161 | 11 | 14.6 |

`tests/test_theorem.py` checks the three signs, the two near-ties, the 8-versus-1 public window, and the 4-versus-1 scaled window.
