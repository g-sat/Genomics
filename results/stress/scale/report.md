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
| heuristic_fixed | 0.301 | 0.541 | 0.711 | 0.604 | 0.399 | 1.004 |
| variant_oracle | 0.276 | 0.686 | 0.700 | 0.684 | 0.310 | 1.001 |
| family_oracle | 0.299 | 0.539 | 0.710 | 0.604 | 0.399 | 1.000 |
| generic_oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 3.624 |
| pedigree_oracle | 0.299 | 0.539 | 0.710 | 0.604 | 0.399 | 1.000 |
| full_oracle | 0.276 | 0.584 | 0.700 | 0.646 | 0.360 | 1.000 |
| random_budget | 0.276 | 0.584 | 0.700 | 0.646 | 0.360 | 1.001 |
| analytical_budget | 0.276 | 0.584 | 0.700 | 0.646 | 0.360 | 1.000 |
| pedpack | 0.276 | 0.584 | 0.700 | 0.646 | 0.360 | 1.000 |
| pedpack_no_pedigree_features | 0.276 | 0.686 | 0.700 | 0.684 | 0.310 | 1.001 |

## Fit-mix execution (median)

| Method | Objective | Mean seconds | Mean bytes | Remote fraction | Mean shuffle |
| --- | ---: | ---: | ---: | ---: | ---: |
| heuristic_fixed | 15.029 | 0.224 | 1.951e+07 | 0.597 | 3.607e+06 |
| variant_oracle | 14.981 | 0.302 | 3.075e+07 | 0.623 | 0.000 |
| family_oracle | 14.964 | 0.223 | 1.951e+07 | 0.609 | 3.607e+06 |
| generic_oracle | 54.228 | 1.053 | 2.470e+07 | 0.621 | 8.100e+07 |
| pedigree_oracle | 14.964 | 0.223 | 1.951e+07 | 0.609 | 3.607e+06 |
| full_oracle | 14.964 | 0.223 | 1.951e+07 | 0.609 | 3.148e+06 |
| random_budget | 14.964 | 0.223 | 1.951e+07 | 0.606 | 3.148e+06 |
| analytical_budget | 14.964 | 0.223 | 1.951e+07 | 0.609 | 3.148e+06 |
| pedpack | 14.964 | 0.223 | 1.951e+07 | 0.609 | 3.148e+06 |
| pedpack_no_pedigree_features | 14.981 | 0.302 | 3.075e+07 | 0.623 | 0.000 |

## Cost model, layouts held out of the pilot

- Neural seconds MAPE: nan
- Analytical seconds MAPE: 0.221
- Neural objective Spearman: nan
- Analytical objective Spearman: 0.990

## Layouts chosen on the fit mix

- Seed 7: 147 samples, 115 trios, 240 layouts. PedPack chose `variant|g=1|w=200|count:16|group_major|id=2` (trio co-location 1.00, 1 files). Generic oracle chose `generic|g=5|w=1600|count:16|group_major|id=86` (trio co-location 0.06, 1 files).
- Seed 11: 102 samples, 70 trios, 240 layouts. PedPack chose `family|g=33|w=3200|count:16|group_major|id=116` (trio co-location 1.00, 3 files). Generic oracle chose `generic|g=4|w=800|count:16|group_major|id=80` (trio co-location 0.06, 1 files).
- Seed 19: 147 samples, 115 trios, 240 layouts. PedPack chose `family|g=46|w=3200|count:16|group_major|id=116` (trio co-location 1.00, 3 files). Generic oracle chose `generic|g=19|w=3200|count:16|group_major|id=44` (trio co-location 0.00, 2 files).

## Amortized cost on the fit mix

Horizon H = 32. Fit shares: {'scan': 0.35, 'family': 0.4, 'extended': 0.25}. Storage weight prices the reference layout at 0.25 reference scans. Amortized seconds per query is the fit objective divided by H. Per sample divides that by the number of people in the cohort.

| Method | Build (s) | Storage (s) | Query (s) | Amortized s/query | Amortized s/sample |
| --- | ---: | ---: | ---: | ---: | ---: |
| pedpack | 3.444 | 0.205 | 0.332 | 0.468 | 3.358e-03 |
| pedpack_no_pedigree_features | 3.444 | 0.205 | 0.354 | 0.468 | 3.420e-03 |
| variant_oracle | 3.444 | 0.205 | 0.354 | 0.468 | 3.420e-03 |
| generic_oracle | 5.211 | 0.205 | 1.522 | 1.695 | 0.012 |
| heuristic_fixed | 4.094 | 0.258 | 0.334 | 0.470 | 3.373e-03 |

## Per-seed fit objectives

| Seed | heuristic_fixed | variant_oracle | family_oracle | generic_oracle | pedigree_oracle | full_oracle | random_budget | analytical_budget | pedpack | pedpack_no_pedigree_features |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 19.408 | 16.480 | 19.344 | 61.585 | 19.344 | 16.480 | 16.499 | 16.480 | 16.480 | 16.480 |
| 11 | 11.008 | 11.162 | 10.961 | 36.610 | 10.961 | 10.961 | 11.151 | 10.961 | 10.961 | 11.162 |
| 19 | 15.029 | 14.981 | 14.964 | 54.228 | 14.964 | 14.964 | 14.964 | 14.964 | 14.964 | 14.981 |

Held-out family and extended queries were not used to choose `pedpack`.
`best_on_heldout_family` in the JSON is an analysis upper bound that did see them.

## Pedigree gap

Same measured fit objective. A gap is `pedpack` choosing a different layout
from `pedpack_no_pedigree_features`, with a lower cost.

- Seed 7: pedpack `variant|g=1|w=200|count:16|group_major|id=2` objective 16.480; no-pedigree `variant|g=1|w=200|count:16|group_major|id=2` objective 16.480; variant `variant|g=1|w=200|count:16|group_major|id=2` objective 16.480; same layout as ablation: True.
- Seed 11: pedpack `family|g=33|w=3200|count:16|group_major|id=116` objective 10.961; no-pedigree `variant|g=1|w=200|count:16|group_major|id=2` objective 11.162; variant `variant|g=1|w=200|count:16|group_major|id=2` objective 11.162; same layout as ablation: False.
- Seed 19: pedpack `family|g=46|w=3200|count:16|group_major|id=116` objective 14.964; no-pedigree `variant|g=1|w=200|count:16|group_major|id=2` objective 14.981; variant `variant|g=1|w=200|count:16|group_major|id=2` objective 14.981; same layout as ablation: False.

## Mixed scan + family evaluation

This list was not used to choose the layout. Its objective uses the same
horizon, with the shares in `mixed_mixture` (default half scan, half family).

- Seed 7: pedpack 19.053, no-pedigree 19.053, generic 79.913.
- Seed 11: pedpack 12.495, no-pedigree 11.452, generic 41.816.
- Seed 19: pedpack 17.651, no-pedigree 15.959, generic 65.616.
