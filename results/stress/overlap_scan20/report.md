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
| heuristic_fixed | 0.359 | 0.547 | 0.686 | 0.622 | 0.448 | 1.033 |
| variant_oracle | 0.347 | 0.654 | 0.691 | 0.681 | 0.294 | 1.000 |
| family_oracle | 0.355 | 0.544 | 0.684 | 0.620 | 0.447 | 1.023 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 2.878 |
| pedigree_oracle | 0.355 | 0.544 | 0.684 | 0.620 | 0.447 | 1.023 |
| full_oracle | 0.347 | 0.654 | 0.680 | 0.666 | 0.294 | 1.000 |
| random_budget | 0.348 | 0.654 | 0.680 | 0.666 | 0.294 | 1.003 |
| analytical_budget | 0.347 | 0.654 | 0.680 | 0.666 | 0.294 | 1.000 |
| pedpack | 0.347 | 0.654 | 0.680 | 0.666 | 0.294 | 1.000 |
| pedpack_no_pedigree_features | 0.347 | 0.654 | 0.691 | 0.681 | 0.294 | 1.000 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 4.552 | 0.103 | 6.359e+06 | 0.566 | 3.480e+06 |
| variant_oracle | 4.406 | 0.099 | 1.005e+07 | 0.619 | 0.000 |
| family_oracle | 4.507 | 0.101 | 6.357e+06 | 0.587 | 3.480e+06 |
| generic_oracle | 12.681 | 0.359 | 1.010e+07 | 0.612 | 2.600e+07 |
| pedigree_oracle | 4.507 | 0.101 | 6.357e+06 | 0.587 | 3.480e+06 |
| full_oracle | 4.406 | 0.099 | 1.005e+07 | 0.619 | 0.000 |
| random_budget | 4.419 | 0.099 | 1.005e+07 | 0.617 | 0.000 |
| analytical_budget | 4.406 | 0.099 | 1.005e+07 | 0.619 | 0.000 |
| pedpack | 4.406 | 0.099 | 1.005e+07 | 0.619 | 0.000 |
| pedpack_no_pedigree_features | 4.406 | 0.099 | 1.005e+07 | 0.619 | 0.000 |

## Cost model, layouts held out of the pilot

- Neural seconds MAPE: nan
- Analytical seconds MAPE: 0.217
- Neural objective Spearman: nan
- Analytical objective Spearman: 0.983

## Layouts chosen on the fit mix

- Seed 7: 104 samples, 84 trios, 240 layouts. PedPack chose `variant|g=1|w=100|count:16|group_major|id=2` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=4|w=400|count:16|group_major|id=80` (trio co-location 0.07, 1 files).
- Seed 11: 73 samples, 53 trios, 240 layouts. PedPack chose `family|g=22|w=1600|count:16|group_major|id=116` (trio co-location 1.00, 2 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.21, 1 files).
- Seed 19: 97 samples, 77 trios, 240 layouts. PedPack chose `variant|g=1|w=100|count:16|group_major|id=2` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=4|w=400|count:16|group_major|id=80` (trio co-location 0.06, 1 files).

## Amortized cost on the fit mix

Horizon H = 32. Fit shares: {'scan': 0.2, 'family': 0.5, 'extended': 0.3}. Storage weight prices the reference layout at 0.25 reference scans. Amortized seconds per query is the fit objective divided by H. Per sample divides that by the number of people in the cohort.

| Method | Build (s) | Storage (s) | Query (s) | Amortized s/query | Amortized s/sample |
| --- | ---: | ---: | ---: | ---: | ---: |
| pedpack | 1.129 | 0.067 | 0.098 | 0.138 | 1.420e-03 |
| pedpack_no_pedigree_features | 1.129 | 0.067 | 0.098 | 0.138 | 1.436e-03 |
| variant_oracle | 1.129 | 0.067 | 0.098 | 0.138 | 1.436e-03 |
| generic_oracle | 1.709 | 0.067 | 0.337 | 0.396 | 3.815e-03 |
| heuristic_fixed | 1.365 | 0.087 | 0.094 | 0.142 | 1.433e-03 |

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 4.552 | 4.406 | 4.507 | 12.681 | 4.507 | 4.406 | 4.419 | 4.406 | 4.406 | 4.406 |
| 11 | 3.348 | 3.502 | 3.316 | 8.911 | 3.316 | 3.316 | 3.417 | 3.316 | 3.316 | 3.502 |
| 19 | 4.672 | 4.459 | 4.623 | 14.174 | 4.623 | 4.459 | 4.459 | 4.459 | 4.459 | 4.459 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `variant|g=1|w=100|count:16|group_major|id=2` objective 4.406; no-pedigree `variant|g=1|w=400|count:16|group_major|id=8` objective 4.406; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 4.406; same layout as ablation: False.
- Seed 11: pedpack `family|g=22|w=1600|count:16|group_major|id=116` objective 3.316; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.502; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.502; same layout as ablation: False.
- Seed 19: pedpack `variant|g=1|w=100|count:16|group_major|id=2` objective 4.459; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 4.459; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 4.459; same layout as ablation: True.

## Mixed scan + family evaluation

This list was not used to choose the layout. Its objective uses the same
horizon, with the shares in `mixed_mixture` (default half scan, half family).

- Seed 7: pedpack 6.064, no-pedigree 6.064, generic 24.642.
- Seed 11: pedpack 5.659, no-pedigree 4.265, generic 15.009.
- Seed 19: pedpack 5.341, no-pedigree 5.341, generic 21.424.
