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
| heuristic_fixed | 0.425 | 0.520 | 0.630 | 0.543 | 0.542 | 1.008 |
| variant_oracle | 0.613 | 0.728 | 0.757 | 0.711 | 0.558 | 1.426 |
| family_oracle | 0.422 | 0.519 | 0.630 | 0.543 | 0.540 | 1.000 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 2.372 |
| pedigree_oracle | 0.422 | 0.519 | 0.630 | 0.543 | 0.540 | 1.000 |
| full_oracle | 0.422 | 0.519 | 0.630 | 0.543 | 0.540 | 1.000 |
| random_budget | 0.422 | 0.519 | 0.630 | 0.543 | 0.540 | 1.000 |
| analytical_budget | 0.422 | 0.519 | 0.630 | 0.543 | 0.540 | 1.000 |
| pedpack | 0.422 | 0.519 | 0.630 | 0.543 | 0.540 | 1.000 |
| pedpack_no_pedigree_features | 0.422 | 0.519 | 0.630 | 0.543 | 0.540 | 1.000 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 600.229 | 8.360 | 8.424e+08 | 0.608 | 0.000 |
| variant_oracle | 889.054 | 18.667 | 1.902e+09 | 0.625 | 0.000 |
| family_oracle | 596.465 | 8.274 | 8.421e+08 | 0.620 | 0.000 |
| generic_oracle | 1.426e+03 | 26.323 | 1.590e+09 | 0.625 | 1.072e+09 |
| pedigree_oracle | 595.726 | 8.260 | 8.419e+08 | 0.618 | 0.000 |
| full_oracle | 595.726 | 8.260 | 8.419e+08 | 0.618 | 0.000 |
| random_budget | 595.726 | 8.260 | 8.419e+08 | 0.618 | 0.000 |
| analytical_budget | 595.726 | 8.260 | 8.419e+08 | 0.618 | 0.000 |
| pedpack | 595.726 | 8.260 | 8.419e+08 | 0.618 | 0.000 |
| pedpack_no_pedigree_features | 595.726 | 8.260 | 8.419e+08 | 0.619 | 0.000 |

## Cost model

- Neural ranker: not run. The selector is the analytical model.
- Analytical seconds MAPE: 0.289
- Analytical objective Spearman: 0.951

## Layouts chosen on the fit mix

- Seed 7: 3375 samples, 636 trios, 108 layouts. PedPack chose `family|g=2110|w=4000|count:32|group_major|id=62` (trio co-location 1.00, 132 files). Generic oracle chose `generic|g=53|w=8000|count:32|group_major|id=32` (trio co-location 0.00, 2 files).
- Seed 11: 3375 samples, 636 trios, 108 layouts. PedPack chose `family|g=2110|w=4000|count:32|group_major|id=62` (trio co-location 1.00, 132 files). Generic oracle chose `generic|g=53|w=8000|count:32|group_major|id=32` (trio co-location 0.00, 2 files).
- Seed 19: 3375 samples, 636 trios, 108 layouts. PedPack chose `pedigree|g=2107|w=8000|count:32|group_major|id=104` (trio co-location 1.00, 66 files). Generic oracle chose `generic|g=53|w=8000|count:32|group_major|id=32` (trio co-location 0.00, 2 files).

## Amortized cost on the fit mix

Horizon H = 32. Fit shares: {'scan': 0.35, 'family': 0.4, 'extended': 0.25}. Storage weight prices the reference layout at 0.25 reference scans. Amortized seconds per query is the fit objective divided by H. Per sample divides that by the number of people in the cohort.

| Method | Build (s) | Storage (s) | Query (s) | Amortized s/query | Amortized s/sample |
| --- | ---: | ---: | ---: | ---: | ---: |
| pedpack | 192.390 | 11.316 | 12.254 | 18.616 | 5.516e-03 |
| pedpack_no_pedigree_features | 192.295 | 11.309 | 12.254 | 18.616 | 5.516e-03 |
| variant_oracle | 199.795 | 11.922 | 21.167 | 27.783 | 8.232e-03 |
| generic_oracle | 307.796 | 11.922 | 34.578 | 44.569 | 0.013 |
| heuristic_fixed | 192.390 | 11.316 | 12.391 | 18.757 | 5.558e-03 |

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 649.397 | 915.094 | 643.882 | 1.527e+03 | 643.882 | 643.882 | 644.718 | 643.882 | 643.882 | 644.853 |
| 11 | 581.169 | 823.070 | 577.231 | 1.342e+03 | 577.231 | 577.231 | 577.269 | 577.231 | 577.231 | 577.269 |
| 19 | 600.229 | 889.054 | 596.465 | 1.426e+03 | 595.726 | 595.726 | 595.726 | 595.726 | 595.726 | 595.726 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `family|g=2110|w=4000|count:32|group_major|id=62` objective 643.882; no-pedigree `family|g=2110|w=8000|count:32|group_major|id=68` objective 644.853; variant `variant|g=1|w=800|count:32|group_major|id=2` objective 915.094; same layout as ablation: False.
- Seed 11: pedpack `family|g=2110|w=4000|count:32|group_major|id=62` objective 577.231; no-pedigree `pedigree|g=2107|w=8000|count:32|group_major|id=104` objective 577.269; variant `variant|g=1|w=800|count:32|group_major|id=2` objective 823.070; same layout as ablation: False.
- Seed 19: pedpack `pedigree|g=2107|w=8000|count:32|group_major|id=104` objective 595.726; no-pedigree `pedigree|g=2107|w=8000|count:32|group_major|id=104` objective 595.726; variant `variant|g=1|w=800|count:32|group_major|id=2` objective 889.054; same layout as ablation: True.
