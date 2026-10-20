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
| heuristic_fixed | 0.340 | 0.603 | 0.724 | 0.655 | 0.470 | 1.243 |
| variant_oracle | 0.274 | 0.703 | 0.711 | 0.693 | 0.329 | 1.000 |
| family_oracle | 0.340 | 0.603 | 0.724 | 0.655 | 0.470 | 1.241 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 3.654 |
| pedigree_oracle | 0.340 | 0.603 | 0.724 | 0.655 | 0.470 | 1.241 |
| full_oracle | 0.274 | 0.703 | 0.711 | 0.693 | 0.329 | 1.000 |
| random_budget | 0.274 | 0.703 | 0.711 | 0.693 | 0.329 | 1.000 |
| analytical_budget | 0.274 | 0.703 | 0.711 | 0.693 | 0.329 | 1.000 |
| pedpack | 0.274 | 0.703 | 0.711 | 0.693 | 0.329 | 1.000 |
| pedpack_no_pedigree_features | 0.274 | 0.703 | 0.711 | 0.693 | 0.329 | 1.000 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 219.012 | 3.533 | 1.942e+08 | 0.614 | 1.558e+08 |
| variant_oracle | 173.136 | 3.501 | 3.567e+08 | 0.625 | 0.000 |
| family_oracle | 218.730 | 3.526 | 1.943e+08 | 0.621 | 1.558e+08 |
| generic_oracle | 642.186 | 12.770 | 3.514e+08 | 0.625 | 9.389e+08 |
| pedigree_oracle | 218.730 | 3.526 | 1.943e+08 | 0.621 | 1.558e+08 |
| full_oracle | 173.136 | 3.501 | 3.567e+08 | 0.625 | 0.000 |
| random_budget | 173.136 | 3.501 | 3.567e+08 | 0.625 | 0.000 |
| analytical_budget | 173.136 | 3.501 | 3.567e+08 | 0.625 | 0.000 |
| pedpack | 173.136 | 3.501 | 3.567e+08 | 0.625 | 0.000 |
| pedpack_no_pedigree_features | 173.136 | 3.501 | 3.567e+08 | 0.625 | 0.000 |

## Cost model, layouts held out of the pilot

- Neural seconds MAPE: nan
- Analytical seconds MAPE: 0.306
- Neural objective Spearman: nan
- Analytical objective Spearman: 0.973

## Layouts chosen on the fit mix

- Seed 7: 667 samples, 507 trios, 108 layouts. PedPack chose `variant|g=1|w=800|count:32|group_major|id=2` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=6|w=4000|count:32|group_major|id=44` (trio co-location 0.02, 1 files).
- Seed 11: 675 samples, 515 trios, 108 layouts. PedPack chose `variant|g=1|w=800|count:32|group_major|id=2` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=6|w=4000|count:32|group_major|id=44` (trio co-location 0.04, 1 files).
- Seed 19: 685 samples, 525 trios, 108 layouts. PedPack chose `variant|g=1|w=800|count:32|group_major|id=2` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=6|w=4000|count:32|group_major|id=44` (trio co-location 0.03, 1 files).

## Amortized cost on the fit mix

Horizon H = 32. Fit shares: {'scan': 0.35, 'family': 0.4, 'extended': 0.25}. Storage weight prices the reference layout at 0.25 reference scans. Amortized seconds per query is the fit objective divided by H. Per sample divides that by the number of people in the cohort.

| Method | Build (s) | Storage (s) | Query (s) | Amortized s/query | Amortized s/sample |
| --- | ---: | ---: | ---: | ---: | ---: |
| pedpack | 39.885 | 2.378 | 4.105 | 5.410 | 8.112e-03 |
| pedpack_no_pedigree_features | 39.885 | 2.378 | 4.105 | 5.410 | 8.112e-03 |
| variant_oracle | 39.885 | 2.378 | 4.105 | 5.410 | 8.112e-03 |
| generic_oracle | 60.565 | 2.379 | 18.111 | 20.068 | 0.030 |
| heuristic_fixed | 47.745 | 3.021 | 5.291 | 6.844 | 0.010 |

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 219.012 | 173.136 | 218.730 | 642.186 | 218.730 | 173.136 | 173.136 | 173.136 | 173.136 | 173.136 |
| 11 | 227.228 | 182.811 | 226.889 | 668.074 | 226.889 | 182.811 | 182.811 | 182.811 | 182.811 | 182.811 |
| 19 | 193.209 | 172.014 | 192.937 | 593.600 | 192.937 | 172.014 | 172.014 | 172.014 | 172.014 | 172.014 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `variant|g=1|w=800|count:32|group_major|id=2` objective 173.136; no-pedigree `variant|g=1|w=800|count:32|group_major|id=2` objective 173.136; variant `variant|g=1|w=800|count:32|group_major|id=2` objective 173.136; same layout as ablation: True.
- Seed 11: pedpack `variant|g=1|w=800|count:32|group_major|id=2` objective 182.811; no-pedigree `variant|g=1|w=800|count:32|group_major|id=2` objective 182.811; variant `variant|g=1|w=800|count:32|group_major|id=2` objective 182.811; same layout as ablation: True.
- Seed 19: pedpack `variant|g=1|w=800|count:32|group_major|id=2` objective 172.014; no-pedigree `variant|g=1|w=800|count:32|group_major|id=2` objective 172.014; variant `variant|g=1|w=800|count:32|group_major|id=2` objective 172.014; same layout as ablation: True.
