# Simulator report

Numbers are from the block simulator in this repository. They are not
HDFS or Hail measurements. A ratio below 1 means a lower objective than
the reference on that row. The fit column amortizes build and storage
over a stated mix of scans and family-scoped queries. Other columns add
the full build to that workload's own query times.

Seeds: 7, 11, 19.

## Objective ratios

Median across seeds. `generic` is the best index-chunk layout after an
exhaustive search on the fit mix. `oracle` is the best layout in the
whole pool on the fit mix, then frozen.

| Method | Fit / generic | Family / generic | Extended / generic | Held-out family / generic | Scan / generic | Fit / oracle |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 0.286 | 0.521 | 0.635 | 0.576 | 0.407 | 1.009 |
| variant_oracle | 0.328 | 0.679 | 0.697 | 0.705 | 0.367 | 1.141 |
| family_oracle | 0.282 | 0.519 | 0.627 | 0.576 | 0.406 | 1.000 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 3.544 |
| pedigree_oracle | 0.282 | 0.519 | 0.627 | 0.576 | 0.406 | 1.000 |
| full_oracle | 0.282 | 0.519 | 0.627 | 0.576 | 0.406 | 1.000 |
| random_budget | 0.287 | 0.524 | 0.627 | 0.573 | 0.414 | 1.001 |
| analytical_budget | 0.282 | 0.519 | 0.627 | 0.576 | 0.406 | 1.000 |
| pedpack | 0.282 | 0.519 | 0.627 | 0.576 | 0.406 | 1.000 |
| pedpack_no_pedigree_features | 0.328 | 0.679 | 0.697 | 0.705 | 0.367 | 1.141 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 3.307 | 0.050 | 5.077e+06 | 0.570 | 0.000 |
| variant_oracle | 3.774 | 0.076 | 7.670e+06 | 0.616 | 0.000 |
| family_oracle | 3.307 | 0.050 | 5.077e+06 | 0.585 | 0.000 |
| generic_oracle | 11.521 | 0.234 | 7.613e+06 | 0.613 | 1.587e+07 |
| pedigree_oracle | 3.307 | 0.050 | 5.077e+06 | 0.585 | 0.000 |
| full_oracle | 3.307 | 0.050 | 5.077e+06 | 0.585 | 0.000 |
| random_budget | 3.309 | 0.050 | 5.101e+06 | 0.584 | 0.000 |
| analytical_budget | 3.307 | 0.050 | 5.077e+06 | 0.585 | 0.000 |
| pedpack | 3.307 | 0.050 | 5.077e+06 | 0.585 | 0.000 |
| pedpack_no_pedigree_features | 3.774 | 0.076 | 7.670e+06 | 0.616 | 0.000 |

## Cost model, layouts held out of the pilot

- Neural seconds MAPE: nan
- Analytical seconds MAPE: 0.207
- Neural objective Spearman: nan
- Analytical objective Spearman: 0.984

## Layouts chosen on the fit mix

- Seed 7: 76 samples, 56 trios, 240 layouts. PedPack chose `family|g=20|w=400|count:16|group_major|id=104` (trio co-location 1.00, 5 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.21, 1 files).
- Seed 11: 58 samples, 38 trios, 240 layouts. PedPack chose `family|g=16|w=1600|count:16|group_major|id=116` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=2|w=400|count:16|group_major|id=80` (trio co-location 0.18, 1 files).
- Seed 19: 67 samples, 47 trios, 240 layouts. PedPack chose `family|g=19|w=1600|count:16|group_major|id=116` (trio co-location 1.00, 2 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.15, 1 files).

## Amortized cost on the fit mix

Horizon H = 32. Fit shares: {'scan': 0.35, 'family': 0.4, 'extended': 0.25}. Storage weight prices the reference layout at 0.25 reference scans. Amortized seconds per query is the fit objective divided by H. Per sample divides that by the number of people in the cohort.

| Method | Build (s) | Storage (s) | Query (s) | Amortized s/query | Amortized s/sample |
| --- | ---: | ---: | ---: | ---: | ---: |
| pedpack | 0.863 | 0.053 | 0.071 | 0.103 | 1.385e-03 |
| pedpack_no_pedigree_features | 0.783 | 0.047 | 0.088 | 0.118 | 1.677e-03 |
| variant_oracle | 0.783 | 0.047 | 0.088 | 0.118 | 1.677e-03 |
| generic_oracle | 1.148 | 0.047 | 0.319 | 0.360 | 4.907e-03 |
| heuristic_fixed | 0.863 | 0.053 | 0.071 | 0.103 | 1.404e-03 |

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 3.307 | 3.774 | 3.307 | 11.521 | 3.307 | 3.307 | 3.309 | 3.307 | 3.307 | 3.774 |
| 11 | 2.606 | 3.113 | 2.570 | 9.107 | 2.570 | 2.570 | 2.668 | 2.570 | 2.570 | 3.113 |
| 19 | 3.469 | 3.805 | 3.438 | 12.738 | 3.438 | 3.438 | 3.438 | 3.438 | 3.438 | 3.805 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `family|g=20|w=400|count:16|group_major|id=104` objective 3.307; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.774; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.774; same layout as ablation: False.
- Seed 11: pedpack `family|g=16|w=1600|count:16|group_major|id=116` objective 2.570; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.113; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.113; same layout as ablation: False.
- Seed 19: pedpack `family|g=19|w=1600|count:16|group_major|id=116` objective 3.438; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.805; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.805; same layout as ablation: False.

## Mixed scan + family evaluation

This list was not used to choose the layout. Its objective uses the same
horizon, with the shares in `mixed_mixture` (default half scan, half family).

- Seed 7: pedpack 4.244, no-pedigree 4.649, generic 15.667.
- Seed 11: pedpack 3.411, no-pedigree 3.768, generic 12.536.
- Seed 19: pedpack 4.047, no-pedigree 4.160, generic 15.212.
