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
| heuristic_fixed | 0.246 | 0.577 | 0.717 | 0.618 | 0.378 | 1.002 |
| variant_oracle | 0.301 | 0.698 | 0.757 | 0.690 | 0.336 | 1.211 |
| family_oracle | 0.246 | 0.577 | 0.717 | 0.618 | 0.378 | 1.000 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.070 |
| pedigree_oracle | 0.246 | 0.577 | 0.717 | 0.618 | 0.378 | 1.000 |
| full_oracle | 0.246 | 0.577 | 0.717 | 0.618 | 0.378 | 1.000 |
| random_budget | 0.249 | 0.578 | 0.717 | 0.618 | 0.380 | 1.004 |
| analytical_budget | 0.246 | 0.577 | 0.717 | 0.618 | 0.378 | 1.000 |
| pedpack | 0.246 | 0.577 | 0.717 | 0.618 | 0.378 | 1.000 |
| pedpack_no_pedigree_features | 0.301 | 0.698 | 0.757 | 0.690 | 0.336 | 1.211 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 721.139 | 10.226 | 1.040e+09 | 0.615 | 0.000 |
| variant_oracle | 895.344 | 18.552 | 1.891e+09 | 0.625 | 0.000 |
| family_oracle | 719.807 | 10.194 | 1.039e+09 | 0.622 | 0.000 |
| generic_oracle | 2.929e+03 | 57.647 | 1.606e+09 | 0.625 | 4.188e+09 |
| pedigree_oracle | 719.807 | 10.194 | 1.039e+09 | 0.622 | 0.000 |
| full_oracle | 719.807 | 10.194 | 1.039e+09 | 0.622 | 0.000 |
| random_budget | 728.182 | 10.371 | 1.041e+09 | 0.616 | 0.000 |
| analytical_budget | 719.807 | 10.194 | 1.039e+09 | 0.622 | 0.000 |
| pedpack | 719.807 | 10.194 | 1.039e+09 | 0.622 | 0.000 |
| pedpack_no_pedigree_features | 895.344 | 18.552 | 1.891e+09 | 0.625 | 0.000 |

## Cost model

- Neural ranker: not run. The selector is the analytical model.
- Analytical seconds MAPE: 0.274
- Analytical objective Spearman: 0.953

## Layouts chosen on the fit mix

- Seed 7: 3370 samples, 2370 trios, 108 layouts. PedPack chose `family|g=994|w=8000|count:32|group_major|id=68` (trio co-location 1.00, 32 files). Generic oracle chose `generic|g=53|w=8000|count:32|group_major|id=32` (trio co-location 0.00, 2 files).
- Seed 11: 3361 samples, 2361 trios, 108 layouts. PedPack chose `family|g=1000|w=8000|count:32|group_major|id=68` (trio co-location 1.00, 32 files). Generic oracle chose `generic|g=53|w=8000|count:32|group_major|id=32` (trio co-location 0.00, 2 files).
- Seed 19: 3328 samples, 2328 trios, 108 layouts. PedPack chose `family|g=989|w=8000|count:32|group_major|id=68` (trio co-location 1.00, 31 files). Generic oracle chose `generic|g=52|w=8000|count:32|group_major|id=32` (trio co-location 0.00, 2 files).

## Amortized cost on the fit mix

Horizon H = 32. Fit shares: {'scan': 0.35, 'family': 0.4, 'extended': 0.25}. Storage weight prices the reference layout at 0.25 reference scans. Amortized seconds per query is the fit objective divided by H. Per sample divides that by the number of people in the cohort.

| Method | Build (s) | Storage (s) | Query (s) | Amortized s/query | Amortized s/sample |
| --- | ---: | ---: | ---: | ---: | ---: |
| pedpack | 220.306 | 13.606 | 15.184 | 22.494 | 6.675e-03 |
| pedpack_no_pedigree_features | 198.923 | 11.869 | 21.375 | 27.980 | 8.303e-03 |
| variant_oracle | 198.923 | 11.869 | 21.375 | 27.980 | 8.303e-03 |
| generic_oracle | 306.337 | 11.869 | 81.567 | 91.541 | 0.027 |
| heuristic_fixed | 220.306 | 13.606 | 15.226 | 22.536 | 6.687e-03 |

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 721.139 | 895.344 | 719.807 | 2.929e+03 | 719.807 | 719.807 | 728.182 | 719.807 | 719.807 | 895.344 |
| 11 | 795.499 | 961.619 | 794.254 | 3.294e+03 | 794.254 | 794.254 | 795.499 | 794.254 | 794.254 | 961.619 |
| 19 | 702.786 | 847.728 | 701.647 | 2.818e+03 | 701.647 | 701.647 | 704.615 | 701.647 | 701.647 | 847.728 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `family|g=994|w=8000|count:32|group_major|id=68` objective 719.807; no-pedigree `variant|g=1|w=800|count:32|group_major|id=2` objective 895.344; variant `variant|g=1|w=800|count:32|group_major|id=2` objective 895.344; same layout as ablation: False.
- Seed 11: pedpack `family|g=1000|w=8000|count:32|group_major|id=68` objective 794.254; no-pedigree `variant|g=1|w=800|count:32|group_major|id=2` objective 961.619; variant `variant|g=1|w=800|count:32|group_major|id=2` objective 961.619; same layout as ablation: False.
- Seed 19: pedpack `family|g=989|w=8000|count:32|group_major|id=68` objective 701.647; no-pedigree `variant|g=1|w=800|count:32|group_major|id=2` objective 847.728; variant `variant|g=1|w=800|count:32|group_major|id=2` objective 847.728; same layout as ablation: False.
