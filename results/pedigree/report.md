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
| heuristic_fixed | 0.286 | 0.530 | 0.637 | 0.597 | 0.392 | 1.007 |
| variant_oracle | 0.335 | 0.703 | 0.680 | 0.720 | 0.349 | 1.168 |
| family_oracle | 0.286 | 0.527 | 0.636 | 0.597 | 0.392 | 1.000 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 3.498 |
| pedigree_oracle | 0.286 | 0.527 | 0.636 | 0.597 | 0.392 | 1.000 |
| full_oracle | 0.286 | 0.527 | 0.636 | 0.597 | 0.392 | 1.000 |
| random_budget | 0.287 | 0.529 | 0.636 | 0.596 | 0.394 | 1.003 |
| analytical_budget | 0.286 | 0.527 | 0.636 | 0.597 | 0.392 | 1.000 |
| pedpack | 0.286 | 0.527 | 0.636 | 0.597 | 0.392 | 1.000 |
| pedpack_no_pedigree_features | 0.335 | 0.703 | 0.680 | 0.720 | 0.349 | 1.168 |
| pedpack_neural | 0.287 | 0.527 | 0.636 | 0.596 | 0.394 | 1.000 |
| pedpack_neural_blind | 0.287 | 0.527 | 0.640 | 0.599 | 0.397 | 1.037 |
| pedpack_generic_pool | 0.335 | 0.703 | 0.680 | 0.720 | 0.349 | 1.168 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 3.224 | 0.048 | 4.799e+06 | 0.558 | 0.000 |
| variant_oracle | 3.705 | 0.073 | 7.431e+06 | 0.616 | 0.000 |
| family_oracle | 3.224 | 0.048 | 4.799e+06 | 0.580 | 0.000 |
| generic_oracle | 11.277 | 0.225 | 7.383e+06 | 0.612 | 1.500e+07 |
| pedigree_oracle | 3.224 | 0.048 | 4.799e+06 | 0.580 | 0.000 |
| full_oracle | 3.224 | 0.048 | 4.799e+06 | 0.580 | 0.000 |
| random_budget | 3.232 | 0.048 | 4.828e+06 | 0.580 | 0.000 |
| analytical_budget | 3.224 | 0.048 | 4.799e+06 | 0.580 | 0.000 |
| pedpack | 3.224 | 0.048 | 4.799e+06 | 0.580 | 0.000 |
| pedpack_no_pedigree_features | 3.705 | 0.073 | 7.431e+06 | 0.616 | 0.000 |
| pedpack_neural | 3.232 | 0.048 | 4.828e+06 | 0.586 | 0.000 |
| pedpack_neural_blind | 3.363 | 0.051 | 4.770e+06 | 0.544 | 0.000 |
| pedpack_generic_pool | 3.705 | 0.073 | 7.431e+06 | 0.616 | 0.000 |

## Cost model, layouts held out of the pilot

- Neural seconds MAPE: 0.888
- Analytical seconds MAPE: 0.215
- Neural objective Spearman: 0.908
- Analytical objective Spearman: 0.987

## Layouts chosen on the fit mix

- Seed 7: 76 samples, 56 trios, 240 layouts. PedPack chose `family|g=20|w=400|count:16|group_major|id=104` (trio co-location 1.00, 5 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.21, 1 files).
- Seed 11: 58 samples, 38 trios, 240 layouts. PedPack chose `family|g=16|w=1600|count:16|group_major|id=116` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=2|w=400|count:16|group_major|id=80` (trio co-location 0.18, 1 files).
- Seed 19: 67 samples, 47 trios, 240 layouts. PedPack chose `family|g=19|w=800|count:16|group_major|id=110` (trio co-location 1.00, 3 files). Generic oracle chose `generic|g=3|w=400|count:16|group_major|id=80` (trio co-location 0.15, 1 files).

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features | pedpack_neural | pedpack_neural_blind | pedpack_generic_pool |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 3.224 | 3.797 | 3.224 | 11.277 | 3.224 | 3.224 | 3.232 | 3.224 | 3.224 | 3.797 | 3.232 | 3.363 | 3.797 |
| 11 | 2.624 | 3.026 | 2.592 | 9.046 | 2.592 | 2.592 | 2.684 | 2.592 | 2.592 | 3.026 | 2.592 | 2.592 | 3.026 |
| 19 | 3.369 | 3.705 | 3.345 | 12.362 | 3.345 | 3.345 | 3.345 | 3.345 | 3.345 | 3.705 | 3.345 | 3.467 | 3.705 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `family|g=20|w=400|count:16|group_major|id=104` objective 3.224; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.797; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.797; same layout as ablation: False.
- Seed 11: pedpack `family|g=16|w=1600|count:16|group_major|id=116` objective 2.592; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.026; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.026; same layout as ablation: False.
- Seed 19: pedpack `family|g=19|w=800|count:16|group_major|id=110` objective 3.345; no-pedigree `variant|g=1|w=100|count:16|group_major|id=2` objective 3.705; variant `variant|g=1|w=100|count:16|group_major|id=2` objective 3.705; same layout as ablation: False.
