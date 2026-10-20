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

| Method | Fit / generic | Scan / generic | Selective / generic | Held-out / generic | Window / generic | Fit / oracle |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 0.491 | 0.561 | 0.677 | 0.588 | 0.775 | 1.069 |
| variant_oracle | 0.409 | 0.405 | 0.752 | 0.728 | 0.748 | 1.000 |
| family_oracle | 0.464 | 0.521 | 0.680 | 0.594 | 0.775 | 1.004 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 2.445 |
| pedigree_oracle | 0.464 | 0.521 | 0.680 | 0.594 | 0.775 | 1.004 |
| full_oracle | 0.409 | 0.405 | 0.752 | 0.728 | 0.748 | 1.000 |
| random_budget | 0.409 | 0.405 | 0.737 | 0.728 | 0.748 | 1.000 |
| analytical_budget | 0.409 | 0.405 | 0.752 | 0.728 | 0.748 | 1.000 |
| pedpack_generic_pool | 0.409 | 0.405 | 0.752 | 0.728 | 0.748 | 1.000 |
| pedpack_no_pedigree_features | 0.409 | 0.405 | 0.752 | 0.728 | 0.748 | 1.000 |
| pedpack | 0.409 | 0.405 | 0.752 | 0.728 | 0.748 | 1.000 |
| pedpack_tuned_on_scan | 0.409 | 0.405 | 0.752 | 0.728 | 0.748 | 1.000 |
| pedpack_tuned_on_selective | 0.515 | 0.563 | 0.672 | 0.578 | 0.775 | 1.310 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 0.530 | 0.014 | 9.222e+05 | 0.000 | 0.000 |
| variant_oracle | 0.460 | 0.011 | 8.883e+05 | 0.163 | 0.000 |
| family_oracle | 0.496 | 0.013 | 9.535e+05 | 0.150 | 0.000 |
| generic_oracle | 1.124 | 0.033 | 1.053e+06 | 0.093 | 1.890e+06 |
| pedigree_oracle | 0.496 | 0.013 | 9.535e+05 | 0.150 | 0.000 |
| full_oracle | 0.460 | 0.011 | 8.883e+05 | 0.163 | 0.000 |
| random_budget | 0.460 | 0.011 | 8.883e+05 | 0.311 | 0.000 |
| analytical_budget | 0.460 | 0.011 | 8.883e+05 | 0.163 | 0.000 |
| pedpack_generic_pool | 0.460 | 0.011 | 8.883e+05 | 0.163 | 0.000 |
| pedpack_no_pedigree_features | 0.460 | 0.011 | 8.883e+05 | 0.163 | 0.000 |
| pedpack | 0.460 | 0.011 | 8.883e+05 | 0.163 | 0.000 |
| pedpack_tuned_on_scan | 0.460 | 0.011 | 8.883e+05 | 0.163 | 0.000 |
| pedpack_tuned_on_selective | 0.649 | 0.018 | 8.847e+05 | 0.376 | 0.000 |

## Cost model, layouts held out of the pilot

- Neural seconds MAPE: 0.383
- Analytical seconds MAPE: 0.241
- Neural objective Spearman: 0.711
- Analytical objective Spearman: 0.982

## Layouts chosen on the fit mix

- Seed 7: 76 samples, 56 trios, 240 layouts. PedPack chose `variant|g=1|w=800|bytes:1048576|group_major|id=16` (trio co-location 1.00, 2 files). Generic oracle chose `generic|g=3|w=100|bytes:1048576|locus_major|id=77` (trio co-location 0.21, 3 files).
- Seed 11: 58 samples, 38 trios, 240 layouts. PedPack chose `variant|g=1|w=100|count:16|group_major|id=2` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=2|w=100|count:16|locus_major|id=75` (trio co-location 0.18, 2 files).
- Seed 19: 67 samples, 47 trios, 240 layouts. PedPack chose `variant|g=1|w=800|naive|group_major|id=12` (trio co-location 1.00, 2 files). Generic oracle chose `generic|g=3|w=100|bytes:1048576|locus_major|id=77` (trio co-location 0.15, 3 files).

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack_generic_pool | pedpack_no_pedigree_features | pedpack | pedpack_tuned_on_scan | pedpack_tuned_on_selective |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 0.551 | 0.460 | 0.521 | 1.124 | 0.521 | 0.460 | 0.460 | 0.460 | 0.460 | 0.460 | 0.460 | 0.460 | 0.537 |
| 11 | 0.456 | 0.426 | 0.428 | 0.914 | 0.428 | 0.426 | 0.426 | 0.426 | 0.426 | 0.426 | 0.426 | 0.426 | 0.823 |
| 19 | 0.530 | 0.496 | 0.496 | 1.260 | 0.496 | 0.496 | 0.496 | 0.496 | 0.496 | 0.496 | 0.496 | 0.496 | 0.649 |

Held-out and window queries were not used to choose `pedpack`.
`best_on_heldout` in the JSON is an analysis upper bound that did see them.

