# Duplication against selectivity

The surface is `family_versus_variant` at 1,200 samples and 8,000 variants,
using the scaled cluster (400 bytes per genotype, 64 KB blocks, replication 3).
A ratio above 1 means the all-samples layout costs more. Black points on
`heatmap.png` are simulator checks where the nuclear layout won; white
points are checks where all-samples won.

Family wins on 32 of 80 formula cells.

## Simulator checks

| Pedigree | Seed | People | Duplication | Fraction | Variant / family | Family wins |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| two_gen | 7 | 157 | 1.318 | 0.05 | 1.120 | True |
| two_gen | 11 | 139 | 1.259 | 0.05 | 1.167 | True |
| two_gen | 7 | 157 | 1.318 | 0.15 | 1.120 | True |
| two_gen | 11 | 139 | 1.259 | 0.15 | 1.126 | True |
| two_gen | 7 | 157 | 1.318 | 0.40 | 1.003 | True |
| two_gen | 11 | 139 | 1.259 | 0.40 | 1.047 | True |
| three_gen | 7 | 147 | 1.449 | 0.05 | 0.892 | False |
| three_gen | 11 | 102 | 1.333 | 0.05 | 1.041 | True |
| three_gen | 7 | 147 | 1.449 | 0.15 | 0.815 | False |
| three_gen | 11 | 102 | 1.333 | 0.15 | 1.004 | True |
| three_gen | 7 | 147 | 1.449 | 0.40 | 0.814 | False |
| three_gen | 11 | 102 | 1.333 | 0.40 | 0.904 | False |
