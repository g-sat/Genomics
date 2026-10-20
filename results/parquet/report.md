# Packed Parquet, one host

Seed 7 of the scaled cohort. Both layouts are Zstandard Parquet with
row groups of 800 variants. The all-samples file
stores every person once. The family side packs nuclear groups
32 to a file (12 files, duplication
1.300). A query reads the covering groups' columns
and the row groups that overlap its variant window. This is the pack
the search chose. It is not one Hail MatrixTable per family.

| Layout | Build s | Bytes |
| --- | ---: | ---: |
| all-samples | 1.69 | 4538050 |
| packed families | 1.28 | 5993701 |

| Query | Groups | Files | Variant s | Variant compressed B | Packed s | Packed compressed B |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| hold-fam-0 | 18 | 10 | 0.060 | 94611 | 0.246 | 94611 |
| hold-fam-1 | 18 | 8 | 0.026 | 97485 | 0.041 | 98871 |
| hold-fam-2 | 18 | 9 | 0.032 | 81391 | 0.036 | 82491 |
| hold-fam-3 | 18 | 10 | 0.030 | 63882 | 0.068 | 64702 |
| scan-0 | 361 | 12 | 0.071 | 2317964 | 0.099 | 3013470 |

Mean held-out seconds: all-samples 0.037, packed 0.098.
Mean held-out compressed bytes: all-samples 84342, packed 85169.
