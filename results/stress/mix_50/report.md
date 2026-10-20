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
| heuristic_fixed | 0.273 | 0.538 | 0.628 | 0.564 | 0.385 | 1.007 |
| variant_oracle | 0.299 | 0.710 | 0.678 | 0.697 | 0.343 | 1.088 |
| family_oracle | 0.273 | 0.536 | 0.626 | 0.564 | 0.385 | 1.000 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 3.661 |
| pedigree_oracle | 0.273 | 0.536 | 0.626 | 0.564 | 0.385 | 1.000 |
| full_oracle | 0.273 | 0.536 | 0.626 | 0.564 | 0.385 | 1.000 |
| random_budget | 0.274 | 0.536 | 0.626 | 0.564 | 0.387 | 1.005 |
| analytical_budget | 0.273 | 0.536 | 0.626 | 0.564 | 0.385 | 1.000 |
| pedpack | 0.273 | 0.536 | 0.626 | 0.564 | 0.385 | 1.000 |
| pedpack_no_pedigree_features | 0.299 | 0.710 | 0.678 | 0.697 | 0.343 | 1.088 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 4.006 | 0.065 | 6.636e+06 | 0.574 | 0.000 |
| variant_oracle | 4.393 | 0.092 | 9.372e+06 | 0.617 | 0.000 |
| family_oracle | 4.006 | 0.065 | 6.636e+06 | 0.590 | 0.000 |
| generic_oracle | 14.665 | 0.309 | 9.294e+06 | 0.615 | 2.162e+07 |
| pedigree_oracle | 4.006 | 0.065 | 6.636e+06 | 0.590 | 0.000 |
| full_oracle | 4.006 | 0.065 | 6.636e+06 | 0.590 | 0.000 |
| random_budget | 4.025 | 0.065 | 6.674e+06 | 0.595 | 0.000 |
| analytical_budget | 4.006 | 0.065 | 6.636e+06 | 0.590 | 0.000 |
| pedpack | 4.006 | 0.065 | 6.636e+06 | 0.590 | 0.000 |
| pedpack_no_pedigree_features | 4.393 | 0.092 | 9.372e+06 | 0.617 | 0.000 |

## Cost model, layouts held out of the pilot

- Neural seconds MAPE: nan
- Analytical seconds MAPE: 0.206
- Neural objective Spearman: nan
- Analytical objective Spearman: 0.989

## Layouts chosen on the fit mix

- Seed 7: 76 samples, 56 trios, 240 layouts. PedPack chose `family|g=20|w=400|count:16|group_major|id=104` (trio co-location 1.00, 5 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.21, 1 files).
- Seed 11: 58 samples, 38 trios, 240 layouts. PedPack chose `family|g=16|w=1600|count:16|group_major|id=116` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=2|w=400|count:16|group_major|id=80` (trio co-location 0.18, 1 files).
- Seed 19: 67 samples, 47 trios, 240 layouts. PedPack chose `family|g=19|w=800|count:16|group_major|id=110` (trio co-location 1.00, 3 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.15, 1 files).

## Amortized cost on the fit mix

Horizon H = 32. Fit shares: {'scan': 0.5, 'family': 0.5, 'extended': 0.0}. Storage weight prices the reference layout at 0.25 reference scans. Amortized seconds per query is the fit objective divided by H. Per sample divides that by the number of people in the cohort.

| Method | Build (s) | Storage (s) | Query (s) | Amortized s/query | Amortized s/sample |
| --- | ---: | ---: | ---: | ---: | ---: |
| pedpack | 0.863 | 0.053 | 0.093 | 0.125 | 1.730e-03 |
| pedpack_no_pedigree_features | 0.783 | 0.047 | 0.108 | 0.137 | 1.882e-03 |
| variant_oracle | 0.783 | 0.047 | 0.108 | 0.137 | 1.882e-03 |
| generic_oracle | 1.148 | 0.047 | 0.417 | 0.458 | 6.296e-03 |
| heuristic_fixed | 0.863 | 0.053 | 0.093 | 0.125 | 1.751e-03 |

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 4.006 | 4.393 | 4.006 | 14.665 | 4.006 | 4.006 | 4.025 | 4.006 | 4.006 | 4.393 |
| 11 | 3.250 | 3.493 | 3.211 | 11.686 | 3.211 | 3.211 | 3.329 | 3.211 | 3.211 | 3.493 |
| 19 | 4.266 | 4.410 | 4.238 | 16.281 | 4.238 | 4.238 | 4.238 | 4.238 | 4.238 | 4.410 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `family|g=20|w=400|count:16|group_major|id=104` objective 4.006; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 4.393; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 4.393; same layout as ablation: False.
- Seed 11: pedpack `family|g=16|w=1600|count:16|group_major|id=116` objective 3.211; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.493; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.493; same layout as ablation: False.
- Seed 19: pedpack `family|g=19|w=800|count:16|group_major|id=110` objective 4.238; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 4.410; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 4.410; same layout as ablation: False.

## Mixed scan + family evaluation

This list was not used to choose the layout. Its objective uses the same
horizon, with the shares in `mixed_mixture` (default half scan, half family).

- Seed 7: pedpack 4.065, no-pedigree 4.499, generic 14.884.
- Seed 11: pedpack 3.010, no-pedigree 3.344, generic 10.803.
- Seed 19: pedpack 3.859, no-pedigree 3.982, generic 14.300.
