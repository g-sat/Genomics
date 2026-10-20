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
| heuristic_fixed | 0.289 | 0.476 | 0.633 | 0.571 | 0.389 | 1.010 |
| variant_oracle | 0.307 | 0.542 | 0.694 | 0.629 | 0.346 | 1.066 |
| family_oracle | 0.288 | 0.470 | 0.630 | 0.570 | 0.389 | 1.000 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 3.474 |
| pedigree_oracle | 0.288 | 0.470 | 0.630 | 0.570 | 0.389 | 1.000 |
| full_oracle | 0.288 | 0.470 | 0.630 | 0.570 | 0.389 | 1.000 |
| random_budget | 0.289 | 0.481 | 0.633 | 0.569 | 0.390 | 1.001 |
| analytical_budget | 0.288 | 0.470 | 0.630 | 0.570 | 0.389 | 1.000 |
| pedpack | 0.288 | 0.470 | 0.630 | 0.570 | 0.389 | 1.000 |
| pedpack_no_pedigree_features | 0.307 | 0.542 | 0.694 | 0.629 | 0.346 | 1.066 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 3.559 | 0.059 | 5.945e+06 | 0.562 | 0.000 |
| variant_oracle | 3.701 | 0.073 | 7.413e+06 | 0.617 | 0.000 |
| family_oracle | 3.549 | 0.059 | 5.959e+06 | 0.572 | 0.000 |
| generic_oracle | 12.306 | 0.261 | 7.514e+06 | 0.611 | 1.850e+07 |
| pedigree_oracle | 3.549 | 0.059 | 5.959e+06 | 0.572 | 0.000 |
| full_oracle | 3.549 | 0.059 | 5.959e+06 | 0.572 | 0.000 |
| random_budget | 3.551 | 0.059 | 5.974e+06 | 0.572 | 0.000 |
| analytical_budget | 3.549 | 0.059 | 5.959e+06 | 0.572 | 0.000 |
| pedpack | 3.549 | 0.059 | 5.959e+06 | 0.572 | 0.000 |
| pedpack_no_pedigree_features | 3.701 | 0.073 | 7.413e+06 | 0.617 | 0.000 |

## Cost model, layouts held out of the pilot

- Neural seconds MAPE: nan
- Analytical seconds MAPE: 0.192
- Neural objective Spearman: nan
- Analytical objective Spearman: 0.988

## Layouts chosen on the fit mix

- Seed 7: 76 samples, 56 trios, 240 layouts. PedPack chose `family|g=20|w=800|count:16|group_major|id=110` (trio co-location 1.00, 3 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.21, 1 files).
- Seed 11: 58 samples, 38 trios, 240 layouts. PedPack chose `family|g=16|w=1600|count:16|group_major|id=116` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=2|w=400|count:16|group_major|id=80` (trio co-location 0.18, 1 files).
- Seed 19: 67 samples, 47 trios, 240 layouts. PedPack chose `family|g=19|w=800|count:16|group_major|id=110` (trio co-location 1.00, 3 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.15, 1 files).

## Amortized cost on the fit mix

Horizon H = 32. Fit shares: {'scan': 0.35, 'family': 0.4, 'extended': 0.25}. Storage weight prices the reference layout at 0.25 reference scans. Amortized seconds per query is the fit objective divided by H. Per sample divides that by the number of people in the cohort.

| Method | Build (s) | Storage (s) | Query (s) | Amortized s/query | Amortized s/sample |
| --- | ---: | ---: | ---: | ---: | ---: |
| pedpack | 0.863 | 0.053 | 0.079 | 0.111 | 1.499e-03 |
| pedpack_no_pedigree_features | 0.783 | 0.047 | 0.089 | 0.116 | 1.694e-03 |
| variant_oracle | 0.783 | 0.047 | 0.089 | 0.116 | 1.694e-03 |
| generic_oracle | 1.148 | 0.047 | 0.343 | 0.385 | 5.210e-03 |
| heuristic_fixed | 0.863 | 0.053 | 0.079 | 0.111 | 1.529e-03 |

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 3.559 | 3.781 | 3.549 | 12.306 | 3.549 | 3.549 | 3.551 | 3.549 | 3.549 | 3.781 |
| 11 | 2.838 | 3.145 | 2.783 | 9.669 | 2.783 | 2.783 | 2.909 | 2.783 | 2.783 | 3.145 |
| 19 | 3.722 | 3.701 | 3.686 | 13.522 | 3.686 | 3.686 | 3.686 | 3.686 | 3.686 | 3.701 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `family|g=20|w=800|count:16|group_major|id=110` objective 3.549; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.781; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.781; same layout as ablation: False.
- Seed 11: pedpack `family|g=16|w=1600|count:16|group_major|id=116` objective 2.783; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.145; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.145; same layout as ablation: False.
- Seed 19: pedpack `family|g=19|w=800|count:16|group_major|id=110` objective 3.686; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.701; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.701; same layout as ablation: False.

## Mixed scan + family evaluation

This list was not used to choose the layout. Its objective uses the same
horizon, with the shares in `mixed_mixture` (default half scan, half family).

- Seed 7: pedpack 4.441, no-pedigree 4.297, generic 15.839.
- Seed 11: pedpack 3.150, no-pedigree 3.209, generic 11.101.
- Seed 19: pedpack 4.449, no-pedigree 4.220, generic 16.713.
