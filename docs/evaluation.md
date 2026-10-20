# Evaluation protocol

The question is whether a cost-driven, pedigree-aware choice of tile width
and packing beats strong layouts once build, storage, and execution are added
up, and whether that gain survives a workload the selector did not fit.

## What this repository measures

| Factor | What is varied | Where |
| --- | --- | --- |
| Layout method | Variant oracle, nuclear family oracle, generic 2D, fixed heuristic, pedigree oracle, full oracle, random budget, pedigree-aware analytical (`pedpack`), blind analytical (`pedpack_no_pedigree_features`), neural and neural-blind as secondary | `study.run_seed` |
| Tile width and pack | Widths, naive vs count vs byte-fill, group-major vs locus-major | `configs/pedigree_gap.json` |
| Workload | Fit mix (scan / family / extended), scan, family-selective, extended, held-out family, held-out extended | `workloads.py` |
| Cluster | Node count 4/8/16, block size at 1/4×, 1×, and 2×, replication 1 and 3, one failed node | sensitivity pass, layouts held fixed |
| Correctness | Mendelian error count on simulated dosages; every query's samples must be covered | enforced in `run_seed` |

The fit mix used for selection is not a raw sum of the query list. Selection
minimizes build plus storage plus a horizon of 32 queries. In
`configs/pedigree_gap.json` that horizon is 35% scans, 40% family-centric
queries, and 25% extended-pedigree queries.

A family-centric query names 5–15% of the nuclear families (at least one) and
a window covering 20–40% of the variant axis. An extended query starts from
one nuclear family, adds the parents of its members, and with probability 1/2
adds a second overlapping family. Held-out queries use different seed families
and freshly drawn windows. An extended held-out query may still name a parent
who also belongs to a training family. None of the held-out lists are in the
selection objective.

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
| `analytical_budget` | Same search as `pedpack` | Allowed |
| `random_budget` | `B` random layouts | Allowed |
| `pedpack` | Pedigree-aware formula picks `B` layouts, then measurement | Allowed |
| `pedpack_no_pedigree_features` | Same budget, blind formula (no group identity, no kinship discount) | Groups exist, but the ranker cannot see who is in them |
| `pedpack_neural` | Secondary. Pilot plus network shortlist | Allowed |
| `pedpack_neural_blind` | Secondary. Pedigree feature columns zeroed | Allowed, but hidden from the network |
| `pedpack_generic_pool` | Secondary. Same neural budget, pool is variant plus index chunks | No |

If `pedpack` and `pedpack_no_pedigree_features` choose the same layout, the
pedigree claim is not supported. The report states that directly. A gap is a
different layout and a lower measured fit cost, and it has to show up on the
held-out family list as well as on the fit mix.

`experiments/run_stress.py` re-fits that comparison at family fractions 0.05,
0.15, and 0.40, at a 50/50 scan/family mix, on a three-generation pedigree,
and on a larger matrix. The locked run in `results/pedigree/` is not
overwritten. The stress summary is `results/stress/summary.md`.

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
| Real multi-node HDFS, not a placement model | NameNode memory, scheduling noise, and disk behavior are not in the simulator. The compose file is one host. This machine has 16 GB of RAM and was not given 4–8 isolated VMs |
| Cluster of 4, 8, 16 DataNodes on separate machines | The sensitivity pass only reassigns replicas in the model |
| Data at tens of GB and above | The scaled simulator is 1,175–1,207 samples and 8,000 variants. That is not tens of GB. Small-file and block-roundup effects still need files that span many blocks for reasons other than the byte scale factor |
| Public genotypes plus pedigree | 1000 Genomes pedigrees, or a cohort you can legally store. `parse_ped` and `load_dosage_npy` are the loaders. Dosages in the current study are simulated on the pedigree |
| Hail or Parquet on a matrix large enough that file-open cost is not the result | `results/hail_columnar/` is one MatrixTable per family and is not the layout. `results/parquet/` is the 12-file pack. Column projection ties the selective bytes. The 16× gap is the opaque-tile model |
| A real codec | Replace `payload_bytes` with measured compressed sizes. The kinship term is a hypothesis about PBWT-style compression, not a measurement of it |
| Failure of one or two DataNodes on the real cluster | The simulator has a failed-node switch. The hardware experiment is still open |

Held-out workloads, the generic 2D baseline, the fixed heuristic, and the
build-plus-storage objective are already in the simulator protocol. A cluster
paper that drops any of them will be arguing a weaker claim than the one in
`docs/method.md`.

## How to read a result

On a family-centric query, a pedigree layout should read fewer bytes than an
all-samples layout, and it should shuffle less than an index-chunk layout
when trios are split. On a pure scan, duplication and file-open cost can
outweigh grouping. PedPack is doing its job when `pedpack` and
`pedpack_no_pedigree_features` choose different layouts and the pedigree-aware
choice has the lower measured cost on the fit mix and on held-out families.

A ratio below 1 against `generic_oracle` is a win on that workload. A ratio
near 1 against `full_oracle` means the budgeted search found the exhaustive
winner. Quote the median ratio and the seed range. Do not quote a single seed
as the result.
