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
| heuristic_fixed | 0.355 | 0.547 | 0.686 | 0.622 | 0.448 | 1.259 |
| variant_oracle | 0.278 | 0.654 | 0.691 | 0.681 | 0.294 | 1.000 |
| family_oracle | 0.352 | 0.544 | 0.684 | 0.620 | 0.447 | 1.249 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 3.599 |
| pedigree_oracle | 0.352 | 0.549 | 0.690 | 0.620 | 0.447 | 1.249 |
| full_oracle | 0.278 | 0.654 | 0.691 | 0.681 | 0.294 | 1.000 |
| random_budget | 0.279 | 0.655 | 0.691 | 0.681 | 0.294 | 1.000 |
| analytical_budget | 0.278 | 0.654 | 0.691 | 0.681 | 0.294 | 1.000 |
| pedpack | 0.278 | 0.654 | 0.691 | 0.681 | 0.294 | 1.000 |
| pedpack_no_pedigree_features | 0.278 | 0.654 | 0.691 | 0.681 | 0.294 | 1.000 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 6.400 | 0.103 | 6.359e+06 | 0.566 | 3.480e+06 |
| variant_oracle | 5.012 | 0.099 | 1.005e+07 | 0.619 | 0.000 |
| family_oracle | 6.344 | 0.101 | 6.357e+06 | 0.587 | 3.480e+06 |
| generic_oracle | 18.035 | 0.359 | 1.010e+07 | 0.612 | 2.600e+07 |
| pedigree_oracle | 6.344 | 0.101 | 7.162e+06 | 0.596 | 3.480e+06 |
| full_oracle | 5.012 | 0.099 | 1.005e+07 | 0.619 | 0.000 |
| random_budget | 5.027 | 0.099 | 1.005e+07 | 0.618 | 0.000 |
| analytical_budget | 5.012 | 0.099 | 1.005e+07 | 0.619 | 0.000 |
| pedpack | 5.012 | 0.099 | 1.005e+07 | 0.619 | 0.000 |
| pedpack_no_pedigree_features | 5.012 | 0.099 | 1.005e+07 | 0.619 | 0.000 |

## Cost model, layouts held out of the pilot

- Neural seconds MAPE: nan
- Analytical seconds MAPE: 0.217
- Neural objective Spearman: nan
- Analytical objective Spearman: 0.983

## Layouts chosen on the fit mix

- Seed 7: 104 samples, 84 trios, 240 layouts. PedPack chose `variant|g=1|w=400|count:16|group_major|id=8` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=4|w=400|count:16|group_major|id=80` (trio co-location 0.07, 1 files).
- Seed 11: 73 samples, 53 trios, 240 layouts. PedPack chose `variant|g=1|w=100|count:16|group_major|id=2` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.21, 1 files).
- Seed 19: 97 samples, 77 trios, 240 layouts. PedPack chose `variant|g=1|w=100|count:16|group_major|id=2` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=4|w=400|count:16|group_major|id=80` (trio co-location 0.06, 1 files).

## Amortized cost on the fit mix

Horizon H = 32. Fit shares: {'scan': 0.35, 'family': 0.4, 'extended': 0.25}. Storage weight prices the reference layout at 0.25 reference scans. Amortized seconds per query is the fit objective divided by H. Per sample divides that by the number of people in the cohort.

| Method | Build (s) | Storage (s) | Query (s) | Amortized s/query | Amortized s/sample |
| --- | ---: | ---: | ---: | ---: | ---: |
| pedpack | 1.129 | 0.067 | 0.116 | 0.157 | 1.706e-03 |
| pedpack_no_pedigree_features | 1.129 | 0.067 | 0.116 | 0.157 | 1.706e-03 |
| variant_oracle | 1.129 | 0.067 | 0.116 | 0.157 | 1.706e-03 |
| generic_oracle | 1.709 | 0.067 | 0.504 | 0.564 | 5.419e-03 |
| heuristic_fixed | 1.365 | 0.087 | 0.152 | 0.200 | 2.038e-03 |

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 6.400 | 5.012 | 6.344 | 18.035 | 6.344 | 5.012 | 5.027 | 5.012 | 5.012 | 5.012 |
| 11 | 4.761 | 4.043 | 4.732 | 12.585 | 4.726 | 4.043 | 4.043 | 4.043 | 4.043 | 4.043 |
| 19 | 6.668 | 5.294 | 6.612 | 20.770 | 6.612 | 5.294 | 5.294 | 5.294 | 5.294 | 5.294 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `variant|g=1|w=400|count:16|group_major|id=8` objective 5.012; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 5.012; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 5.012; same layout as ablation: False.
- Seed 11: pedpack `variant|g=1|w=100|count:16|group_major|id=2` objective 4.043; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 4.043; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 4.043; same layout as ablation: True.
- Seed 19: pedpack `variant|g=1|w=100|count:16|group_major|id=2` objective 5.294; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 5.294; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 5.294; same layout as ablation: True.

## Mixed scan + family evaluation

This list was not used to choose the layout. Its objective uses the same
horizon, with the shares in `mixed_mixture` (default half scan, half family).

- Seed 7: pedpack 6.064, no-pedigree 6.064, generic 24.642.
- Seed 11: pedpack 4.265, no-pedigree 4.265, generic 15.009.
- Seed 19: pedpack 5.341, no-pedigree 5.341, generic 21.424.
