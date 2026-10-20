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
| heuristic_fixed | 0.286 | 0.539 | 0.638 | 0.603 | 0.392 | 1.007 |
| variant_oracle | 0.339 | 0.741 | 0.692 | 0.747 | 0.353 | 1.185 |
| family_oracle | 0.286 | 0.540 | 0.633 | 0.603 | 0.392 | 1.000 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 3.498 |
| pedigree_oracle | 0.286 | 0.540 | 0.633 | 0.603 | 0.392 | 1.000 |
| full_oracle | 0.286 | 0.540 | 0.633 | 0.603 | 0.392 | 1.000 |
| random_budget | 0.286 | 0.540 | 0.640 | 0.603 | 0.394 | 1.001 |
| analytical_budget | 0.286 | 0.540 | 0.633 | 0.603 | 0.392 | 1.000 |
| pedpack | 0.286 | 0.540 | 0.633 | 0.603 | 0.392 | 1.000 |
| pedpack_no_pedigree_features | 0.339 | 0.741 | 0.692 | 0.747 | 0.353 | 1.185 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 3.193 | 0.047 | 4.729e+06 | 0.564 | 0.000 |
| variant_oracle | 3.712 | 0.074 | 7.460e+06 | 0.616 | 0.000 |
| family_oracle | 3.193 | 0.047 | 4.729e+06 | 0.582 | 0.000 |
| generic_oracle | 11.170 | 0.222 | 7.438e+06 | 0.609 | 1.463e+07 |
| pedigree_oracle | 3.193 | 0.047 | 4.729e+06 | 0.582 | 0.000 |
| full_oracle | 3.193 | 0.047 | 4.729e+06 | 0.582 | 0.000 |
| random_budget | 3.197 | 0.046 | 4.730e+06 | 0.582 | 0.000 |
| analytical_budget | 3.193 | 0.047 | 4.729e+06 | 0.582 | 0.000 |
| pedpack | 3.193 | 0.047 | 4.729e+06 | 0.582 | 0.000 |
| pedpack_no_pedigree_features | 3.712 | 0.074 | 7.460e+06 | 0.616 | 0.000 |

## Cost model, layouts held out of the pilot

- Neural seconds MAPE: nan
- Analytical seconds MAPE: 0.227
- Neural objective Spearman: nan
- Analytical objective Spearman: 0.985

## Layouts chosen on the fit mix

- Seed 7: 76 samples, 56 trios, 240 layouts. PedPack chose `family|g=20|w=400|count:16|group_major|id=104` (trio co-location 1.00, 5 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.21, 1 files).
- Seed 11: 58 samples, 38 trios, 240 layouts. PedPack chose `family|g=16|w=1600|count:16|group_major|id=116` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=2|w=400|count:16|group_major|id=80` (trio co-location 0.18, 1 files).
- Seed 19: 67 samples, 47 trios, 240 layouts. PedPack chose `family|g=19|w=800|count:16|group_major|id=110` (trio co-location 1.00, 3 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.15, 1 files).

## Amortized cost on the fit mix

Horizon H = 32. Fit shares: {'scan': 0.35, 'family': 0.4, 'extended': 0.25}. Storage weight prices the reference layout at 0.25 reference scans. Amortized seconds per query is the fit objective divided by H. Per sample divides that by the number of people in the cohort.

| Method | Build (s) | Storage (s) | Query (s) | Amortized s/query | Amortized s/sample |
| --- | ---: | ---: | ---: | ---: | ---: |
| pedpack | 0.863 | 0.053 | 0.068 | 0.100 | 1.369e-03 |
| pedpack_no_pedigree_features | 0.783 | 0.047 | 0.089 | 0.116 | 1.633e-03 |
| variant_oracle | 0.783 | 0.047 | 0.089 | 0.116 | 1.633e-03 |
| generic_oracle | 1.148 | 0.047 | 0.308 | 0.349 | 4.789e-03 |
| heuristic_fixed | 0.863 | 0.053 | 0.068 | 0.100 | 1.385e-03 |

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 3.193 | 3.784 | 3.193 | 11.170 | 3.193 | 3.193 | 3.197 | 3.193 | 3.193 | 3.784 |
| 11 | 2.571 | 3.031 | 2.541 | 8.889 | 2.541 | 2.541 | 2.637 | 2.541 | 2.541 | 3.031 |
| 19 | 3.335 | 3.712 | 3.312 | 12.242 | 3.312 | 3.312 | 3.312 | 3.312 | 3.312 | 3.712 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `family|g=20|w=400|count:16|group_major|id=104` objective 3.193; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.784; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.784; same layout as ablation: False.
- Seed 11: pedpack `family|g=16|w=1600|count:16|group_major|id=116` objective 2.541; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.031; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.031; same layout as ablation: False.
- Seed 19: pedpack `family|g=19|w=800|count:16|group_major|id=110` objective 3.312; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.712; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.712; same layout as ablation: False.

## Mixed scan + family evaluation

This list was not used to choose the layout. Its objective uses the same
horizon, with the shares in `mixed_mixture` (default half scan, half family).

- Seed 7: pedpack 4.118, no-pedigree 4.588, generic 15.057.
- Seed 11: pedpack 3.281, no-pedigree 3.621, generic 12.150.
- Seed 19: pedpack 3.884, no-pedigree 4.063, generic 14.455.
