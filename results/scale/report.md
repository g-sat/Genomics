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
| heuristic_fixed | 0.244 | 0.564 | 0.702 | 0.612 | 0.378 | 1.002 |
| variant_oracle | 0.297 | 0.704 | 0.729 | 0.689 | 0.338 | 1.218 |
| family_oracle | 0.244 | 0.563 | 0.702 | 0.612 | 0.378 | 1.000 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.097 |
| pedigree_oracle | 0.244 | 0.563 | 0.702 | 0.612 | 0.378 | 1.000 |
| full_oracle | 0.244 | 0.563 | 0.702 | 0.612 | 0.378 | 1.000 |
| random_budget | 0.247 | 0.564 | 0.702 | 0.613 | 0.380 | 1.010 |
| analytical_budget | 0.244 | 0.563 | 0.702 | 0.612 | 0.378 | 1.000 |
| pedpack | 0.244 | 0.563 | 0.702 | 0.612 | 0.378 | 1.000 |
| pedpack_no_pedigree_features | 0.297 | 0.704 | 0.729 | 0.689 | 0.338 | 1.218 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 267.615 | 3.855 | 3.921e+08 | 0.615 | 0.000 |
| variant_oracle | 323.605 | 6.622 | 6.748e+08 | 0.625 | 0.000 |
| family_oracle | 267.199 | 3.845 | 3.920e+08 | 0.622 | 0.000 |
| generic_oracle | 1.095e+03 | 21.798 | 6.109e+08 | 0.625 | 1.562e+09 |
| pedigree_oracle | 267.199 | 3.845 | 3.920e+08 | 0.622 | 0.000 |
| full_oracle | 267.199 | 3.845 | 3.920e+08 | 0.622 | 0.000 |
| random_budget | 269.870 | 3.901 | 3.926e+08 | 0.611 | 0.000 |
| analytical_budget | 267.199 | 3.845 | 3.920e+08 | 0.622 | 0.000 |
| pedpack | 267.199 | 3.845 | 3.920e+08 | 0.622 | 0.000 |
| pedpack_no_pedigree_features | 323.605 | 6.622 | 6.748e+08 | 0.625 | 0.000 |

## Cost model, layouts held out of the pilot

- Neural seconds MAPE: nan
- Analytical seconds MAPE: 0.280
- Neural objective Spearman: nan
- Analytical objective Spearman: 0.948

## Layouts chosen on the fit mix

- Seed 7: 1207 samples, 847 trios, 108 layouts. PedPack chose `family|g=361|w=8000|count:32|group_major|id=68` (trio co-location 1.00, 12 files). Generic oracle chose `generic|g=19|w=8000|count:32|group_major|id=32` (trio co-location 0.00, 1 files).
- Seed 11: 1175 samples, 815 trios, 108 layouts. PedPack chose `family|g=351|w=8000|count:32|group_major|id=68` (trio co-location 1.00, 11 files). Generic oracle chose `generic|g=19|w=8000|count:32|group_major|id=32` (trio co-location 0.00, 1 files).
- Seed 19: 1176 samples, 816 trios, 108 layouts. PedPack chose `family|g=345|w=4000|count:32|group_major|id=62` (trio co-location 1.00, 22 files). Generic oracle chose `generic|g=19|w=8000|count:32|group_major|id=32` (trio co-location 0.00, 1 files).

## Amortized cost on the fit mix

Horizon H = 32. Fit shares: {'scan': 0.35, 'family': 0.4, 'extended': 0.25}. Storage weight prices the reference layout at 0.25 reference scans. Amortized seconds per query is the fit objective divided by H. Per sample divides that by the number of people in the cohort.

| Method | Build (s) | Storage (s) | Query (s) | Amortized s/query | Amortized s/sample |
| --- | ---: | ---: | ---: | ---: | ---: |
| pedpack | 76.773 | 4.741 | 5.720 | 8.350 | 6.918e-03 |
| pedpack_no_pedigree_features | 69.574 | 4.151 | 7.809 | 10.113 | 8.429e-03 |
| variant_oracle | 69.574 | 4.151 | 7.809 | 10.113 | 8.429e-03 |
| generic_oracle | 107.156 | 4.150 | 30.640 | 34.206 | 0.028 |
| heuristic_fixed | 76.773 | 4.741 | 5.733 | 8.363 | 6.929e-03 |

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 267.615 | 325.563 | 267.199 | 1.095e+03 | 267.199 | 267.199 | 269.870 | 267.199 | 267.199 | 325.563 |
| 11 | 270.562 | 323.605 | 270.101 | 1.108e+03 | 270.101 | 270.101 | 270.562 | 270.101 | 270.101 | 323.605 |
| 19 | 220.361 | 282.578 | 220.136 | 873.600 | 220.136 | 220.136 | 223.371 | 220.136 | 220.136 | 282.578 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `family|g=361|w=8000|count:32|group_major|id=68` objective 267.199; no-pedigree `variant|g=1|w=800|count:32|group_major|id=2` objective 325.563; variant `variant|g=1|w=800|count:32|group_major|id=2` objective 325.563; same layout as ablation: False.
- Seed 11: pedpack `family|g=351|w=8000|count:32|group_major|id=68` objective 270.101; no-pedigree `variant|g=1|w=800|count:32|group_major|id=2` objective 323.605; variant `variant|g=1|w=800|count:32|group_major|id=2` objective 323.605; same layout as ablation: False.
- Seed 19: pedpack `family|g=345|w=4000|count:32|group_major|id=62` objective 220.136; no-pedigree `variant|g=1|w=4000|count:32|group_major|id=8` objective 282.578; variant `variant|g=1|w=4000|count:32|group_major|id=8` objective 282.578; same layout as ablation: False.
