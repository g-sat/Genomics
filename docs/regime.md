# Locked regime: the pedigree gap

This is the simulator result in which pedigree information changes the layout.
`results/pedigree/study.json` is that run. Do not overwrite it with a stress
test. Stress tests live in `results/stress/` and start from the same config.

## What was run

| Item | Value |
| --- | --- |
| Config | `configs/pedigree_gap.json` |
| Seeds | 7, 11, 19 |
| Cohort | 10 founder couples, 2 generations, mean 2.2 children, cross rate 0.25, at least 48 people |
| Matrix | 1,600 variants, 100 bp spacing, dosages transmitted without Mendelian error |
| Cluster | 8 nodes, 64 KB blocks, replication 3, 400 bytes per genotype, compression base 0.45, kinship alpha 0.55, file open 0.001 s |
| Horizon | 32 queries |
| Fit shares | scan 0.35, family 0.40, extended 0.25 |
| Storage price | replicated bytes of the reference all-samples layout cost 0.25 reference scans |
| Budget | 36 measured layouts out of the pool |
| Selector | `pedpack` ranks with the pedigree-aware formula, then keeps the cheapest measured layout |
| Ablation | `pedpack_no_pedigree_features` uses the same budget and a formula that sees group sizes only |

The config that produced the JSON is stored again inside that JSON under `config`.

## Queries

Family-centric queries name 5–15% of the nuclear families (at least one) and
a window covering 20–40% of the variant axis. Extended queries add the parents
of a nuclear family and, half the time, a second family. Held-out queries use
different seed families and newly drawn windows. Counts: 4 fit scans, 8 fit
family queries, 6 fit extended queries, and the evaluation lists in the config.

## Layout pool

Widths 100, 400, 800, and the full 1,600 variants. Orders `group_major` and
`locus_major`. Packs: one tile per file, 16 tiles per file, and a byte target
of 4 blocks. Groups: all samples, index chunks of 8, 16, and 32 after column
scrambling, nuclear families, kinship agglomeration at sizes 12 and 24 with
minimum kinship 0.0625 and encoder weights 0 and 0.5, three-generation groups,
and one fixed nuclear heuristic at width 400 and `count:16`.

## Cost model

`pedpack.cost.analytical_query` with `pedigree_aware=True` opens the true
covering groups, uses kinship-compressed tile sizes, rounds each file up to
a block, and charges shuffle from the query's own trios. The remote fraction
is a function of replication and of how many files are opened. It does not
place the task.

The ablation sets `pedigree_aware=False`. A group of size `s` is opened with
probability `1 - (1 - q/N)^s`. Tiles are sized with kinship 0. Co-location is
mean group size over N.

The returned objective is

```
C = T_build + ω S + H (0.35 mean T_scan + 0.40 mean T_family + 0.25 mean T_extended)
```

`C / H` is the amortized seconds of one query from that mix. `(C / H) / N`
is the same quantity per sample. Build and storage are not charged again on
top of that. A column that adds one full build to a short query list is a
different number; the report labels it as such.

## What this run showed

On every seed, PedPack selected a nuclear-family layout and matched the
exhaustive winner. The ablation selected the all-samples layout that the
variant oracle selected. Fit objectives were 3.224 vs 3.797, 2.592 vs 3.026,
and 3.345 vs 3.705. The formula's objective Spearman against the simulator,
on layouts outside the neural pilot, was 0.987. Seconds MAPE was 0.215.

The three re-fits that keep this gap are frozen as
`configs/positive/select_05.json`, `select_15.json`, and `mix_50.json`.
The map of where the gap closes is `docs/boundary.md` and
`results/stress/boundary.png`.

## What is not locked

Stress tests that start from this config and change one factor are in
`results/stress/summary.md`. Selectivity of 5%, 15%, and 40%, and a 50/50
scan/family mix, keep the gap on seeds 7, 11, and 19. Three-generation
pedigrees with duplication around 1.4 do not: both rankers then select the
all-samples layout, and that layout is the exhaustive winner. The same
split, at 1,175–1,207 people and at 667–685 people, is in
`results/scale/` and `results/scale_overlap/`. A real multi-node HDFS
measurement is not this result. The neural ranker is not the method.
