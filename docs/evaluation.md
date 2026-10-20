# Evaluation protocol

The question is whether a cost-driven, pedigree-aware choice of tile width
and packing beats strong layouts once build, storage, and execution are added
up, and whether that gain survives a workload the selector did not fit.

## What this repository measures

| Factor | What is varied | Where |
| --- | --- | --- |
| Layout method | Variant, nuclear family, generic 2D, fixed heuristic, pedigree oracle, full oracle, random budget, analytical budget, PedPack, pedigree features hidden, pedigree groups removed, retuned on scans, retuned on selective queries | `study.run_seed` |
| Tile width and pack | Widths, naive vs count vs byte-fill, group-major vs locus-major | `configs/small.json` |
| Workload | Fit mix, scan, selective, held-out families, short windows | `workloads.py` |
| Cluster | Node count 4/8/16, block size at 1/4×, 1×, and 2×, replication 1 and 3, one failed node | sensitivity pass, layouts held fixed |
| Correctness | Mendelian error count on simulated dosages; every query's samples must be covered | enforced in `run_seed` |

The fit mix used for selection is not a raw sum of the query list. Full scans
read so much more than a family-scoped query that an unweighted sum is a scan
benchmark. Selection minimizes build plus storage plus a horizon of 32 queries,
with 25% of that horizon on the fit scans and 75% on the fit selective queries.
`pedpack_tuned_on_scan` and `pedpack_tuned_on_selective` repeat the search with
those workloads as the objective. Held-out queries use nuclear families that
were withheld from the fit list, including narrow regions and wider ones.
Window queries are short spans on those held-out families. They are an access
pattern, not a variant caller. Neither set is in the selection objective.

Repetition is across pedigree seeds, not across identical reruns of a
deterministic simulator. The report gives the median and the min/max of the
objective. Ratios are computed inside a seed, against that seed's generic
oracle and full oracle, then the median ratio is reported.

## Baselines

| Name | Search | Pedigree groups |
| --- | --- | --- |
| `heuristic_fixed` | None. Nuclear families, one width, one pack | Yes |
| `variant_oracle` | Exhaustive inside the all-samples grouping | No |
| `family_oracle` | Exhaustive inside nuclear families | Yes |
| `generic_oracle` | Exhaustive inside index chunks | No |
| `pedigree_oracle` | Exhaustive inside pedigree groupings | Yes |
| `full_oracle` | Exhaustive on the whole pool | Allowed |
| `analytical_budget` | Formula picks `B` layouts, then measurement | Allowed |
| `random_budget` | `B` random layouts | Allowed |
| `pedpack` | Pilot plus network shortlist, budget `B` | Allowed |
| `pedpack_no_pedigree_features` | Same budget and same pilot, pedigree columns zeroed | Allowed, but hidden from the ranker |
| `pedpack_generic_pool` | Same budget, pool is variant plus index chunks | No |

If `pedpack` does not beat `pedpack_generic_pool` and `generic_oracle` on the
mixed and held-out workloads, the pedigree-aware claim is not supported by
this run. The report is expected to show that, not to hide it. Full scans can
favor a coarse variant layout; `pedpack_tuned_on_scan` exists so that shift is
visible.

## Metrics

For every method and workload the JSON stores objective, mean seconds, the
standard deviation across queries, mean bytes, byte-weighted remote fraction,
mean shuffle, files touched, file count, block count, replicated bytes, build
seconds, duplication, and trio co-location.

The cost network is scored on layouts that were not in its pilot: seconds
MAPE against the simulator, and Spearman correlation of the predicted
objective against the measured objective. The analytical formula is scored on
the same layouts.

## What is not measured here

These are required before the paper can claim a distributed-systems result.
They are not optional polish.

| Gap | Why it matters |
| --- | --- |
| Real multi-node HDFS, not a placement model | NameNode memory, scheduling noise, and disk behavior are not in the simulator |
| Cluster of 4, 8, 16 DataNodes on separate machines | The sensitivity pass only reassigns replicas in the model |
| Data at tens of GB and above | Small-file and block-roundup effects need files that span many blocks for reasons other than the byte scale factor |
| Public genotypes plus pedigree | 1000 Genomes pedigrees, or a cohort you can legally store. `parse_ped` and `load_dosage_npy` are the loaders. Dosages in the current study are simulated on the pedigree |
| Hail on the same queries | MatrixTable partitioned by locus, and a family-keyed table if you build one. Report runtime, bytes read, and shuffle from the Hail log. Do not cite a number this repository did not produce |
| A real codec | Replace `payload_bytes` with measured compressed sizes. The kinship term is a hypothesis about PBWT-style compression, not a measurement of it |
| Failure of one or two DataNodes on the real cluster | The simulator has a failed-node switch. The hardware experiment is still open |

Held-out workloads, the generic 2D baseline, the fixed heuristic, and the
build-plus-storage objective are already in the simulator protocol. A cluster
paper that drops any of them will be arguing a weaker claim than the one in
`docs/method.md`.

## How to read a result

On a mixed fit workload, a pedigree layout should read fewer bytes than an
all-samples layout for selective queries, and it should shuffle less than an
index-chunk layout when trios are split. On a pure scan, file-open cost and
sequential layout can outweigh grouping. PedPack is doing its job when the
chosen layout changes between `pedpack` and `pedpack_tuned_on_scan`.

A ratio below 1 against `generic_oracle` is a win on that workload. A ratio
near 1 against `full_oracle` means the budgeted search found the exhaustive
winner. Quote the median ratio and the seed range. Do not quote a single seed
as the result.
