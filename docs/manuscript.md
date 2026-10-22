# Pedigree-Aware Physical Design for Distributed Genomic Storage

## Abstract

PedPack is a storage policy for an opaque block engine: HDFS, or any blob
store that returns a byte range and cannot project a named sample out of a
tile. It chooses tile width, sample grouping, and file packing for a genotype
matrix whose queries mix full scans with family-scoped hereditary lookups.
The objective adds build time, replicated storage, and expected query time.
The selector ranks candidate layouts with a closed-form cost model that can
see which samples share a group, then keeps the cheapest layout a block
simulator measures inside a fixed budget. A second ranker sees the same group
sizes and the same packing, and does not see who is related to whom. Hail
and Parquet are the check that this policy does not transfer to a columnar
engine, which can already return the named samples from one file.

On 2-generation pedigrees the nuclear layout stores each overlapping person
about 1.21–1.30 times. From 58 people through about 3,370, with family queries
naming about 5–15% of families, PedPack selects that layout and matches an
exhaustive search of the pool. The ranker that cannot see pedigree identity
selects the all-samples layout. At 1,175–3,370 people and 8,000 variants its
fit objective is 1.20–1.28× the family layout, and it reads 16–18× more bytes
on held-out families. The family layout costs more to build and to store.
The win is query time. At nuclear duplication 1.42–1.44 both rankers select
all-samples, and the family oracle is the more expensive layout. A two-layout
formula with the same block model predicts that sign. A neural cost model is
not part of the method. These numbers are from a block simulator. A separate
local-mode Hail run is reported below. Neither is a multi-node HDFS result.

## 1. The design problem

A dosage matrix has M variants and N samples. Hereditary checks need the
child and both parents at the same loci. A scan needs every sample. The
physical design decides, before the queries arrive, how that matrix is cut
and packed.

A logical tile is one sample group over a contiguous run of variants. Tiles
are concatenated into files. A query reads whole tiles, expanded to block
boundaries. If a sample sits in two groups, both copies are stored and the
planner reads one covering set of groups.

Five groupings are in the pool.

- All samples in one group. Every trio is co-located because everyone is.
- Index-contiguous chunks after the columns have been scrambled, so the chunks
  are not secretly families.
- Nuclear families: the two parents and their children. Across two generations
  a person can sit in the family they were born into and the family where they
  are a parent. The copies are intentional and they are billed.
- Kinship agglomeration, size-capped, optionally biased by a kinship encoder,
  plus a three-generation expansion that overlaps.
- A fixed heuristic: nuclear families, width 400 variants, 16 tiles per file.
  No search.

Packing is one tile per file, a fixed tile count per file, or a byte target
of four blocks. Write order is group-major or locus-major.

A trio is co-located when one group contains the father, the mother, and the
child. The layout rate is that fraction over the whole pedigree. The query
rate is that fraction over the trios named in the query.

The reason to search this space, rather than to declare that related samples
belong together, is that "together" has two cheap implementations. One tile
of everyone is together. One tile per nuclear family is together for that
family and expensive for a scan. Which one wins depends on the fraction of
scans, the duplication, the record size, and how many families a query
actually names. The selector exists to make that call from a cost.

## 2. Cost

For a layout L and a horizon of H queries,

```
C(L) = T_build(L) + ω S(L)
     + H (α_scan mean T_scan + α_family mean T_family + α_extended mean T_extended)
```

S(L) is the replicated size. ω is set so that storing the reference
all-samples layout costs 0.25 reference scans. In the reported run
α_scan = 0.35, α_family = 0.40, α_extended = 0.25, and H = 32. Build,
storage, and query time are all in seconds.

An unweighted sum of a query list is the wrong objective here. A scan reads
the matrix; a family-scoped query reads a few percent of the samples. The
mixture makes the scan share a parameter instead of an accident of list
length.

Tile bytes are the raw record size times a compression factor

```
clip(0.45 × (1 − 0.55 × mean within-group kinship), 0.08, 1)
```

That factor is a stand-in for a codec that shrinks on related haplotypes. It
is not a measured PBWT. Build time reads the source once, shuffles a fraction
(1 − trio co-location) of it, writes the primary copy, and writes extra
replicas at remote bandwidth.

The pedigree-aware formula opens the true covering groups, uses those
compressed tile sizes, rounds each touched file up to a block, and charges
shuffle from the trios in the query. The remote fraction is not task
placement. It starts from the replication shortfall and rises with the number
of files opened. The blind formula keeps group sizes and packing. A group of
size s is opened with probability 1 − (1 − q/N)^s, tiles are sized with
kinship 0, and co-location is the mean group size over N. Two layouts with
the same size histogram get the same selective-byte estimate.

The formula still does not place the task. On the layouts PedPack kept, the
simulator's remote fraction is about 0.56–0.59, and the formula's absolute
objective is low by about a quarter to a third (predicted/measured 0.67–0.77
on the chosen layout). Rank correlation against the simulator is what the
search uses. On the scaled pool it is 0.95.

The same block arithmetic, reduced to two layouts, is
`family_versus_variant`. A query that names a fraction f of the nuclear
groups reads about f of the family layout and the whole all-samples layout.
Extended queries count as twice that fraction. Both layouts are charged the
measured remote fraction of 0.6, because the simulator does not make a
one-file layout local. At 1,200 samples, duplication 1.30, and f = 0.05 the
all-samples objective is 1.06× the family objective. At duplication 1.43 and
the same f, the family layout loses (ratio 0.97). At f = 0.60 and duplication
1.30 it loses as well. The 40% family-query cell is the thin edge: the
formula already prefers all-samples, while the small-cohort simulator still
picks families with a fit ratio near 1. The surface is
`results/boundary/heatmap.png`.

## 3. What searches

PedPack ranks the pool with the pedigree-aware formula and measures the top
B layouts (36 in the locked run, 24 in the scaled run). The returned layout
is the cheapest of those measurements. `pedpack_no_pedigree_features` is the
same budget with the blind formula. `analytical_budget` is the same search
as PedPack. `random_budget` measures B random layouts.

Oracles measure every layout in one family of candidates, or the whole pool.
They are analysis bounds, not a system a user can run without paying for the
full measurement. On the selective regimes PedPack's fit objective matches
the full oracle (ratio 1.000): the formula put the winner inside the budget.
On the high-duplication regimes the family oracle is more expensive than
all-samples, and PedPack matches the variant oracle instead. Matching an
oracle there means the search refused the costly layout.

A kinship encoder and a cost network exist in the code. They are not the
selector. Reports that do not train the network no longer print a missing
error as if the network had been scored.

## 4. Protocol

Cohorts are simulated: 10 founder couples, two generations of children, a
25% chance that a couple spans birth families. Dosages are transmitted
without Mendelian error; the study aborts if that check fails. Columns are
permuted so that consecutive samples are not a family. Each seed draws its
own pedigree and its own queries. Seeds are 7, 11, and 19.

A family-centric query names 5–15% of the nuclear families and a window
covering 20–40% of the variant axis. That fraction is a clinical or
rare-disease request: the affected family, the parents, and sometimes a
second family that shares a parent, against a variant window rather than a
genome-wide scan of the biobank. In a cohort of tens to a few hundred
families, one to a few pedigrees is a few percent to about 15% of the
families. A QC or allele-frequency pass is the scan share. The boundary is
the point of the paper. A GWAS-shaped mix is supposed to keep the
all-samples layout. An extended query starts from one nuclear family, adds
the parents of its members, and with probability one half adds a second
family, so the request is not a clean independent trio.
Held-out queries use different seed families and freshly drawn windows. An
extended held-out query may still name a parent who also sits in a training
family. Held-out lists are not in the selection objective.

The simulator places replicas on 8 nodes at replication 3 and a 64 KB block,
with 400 bytes per genotype and a 1 ms file-open cost. It reads a block from
the node that holds the most replicas of the file's touched blocks, and
charges shuffle when a queried trio is split across the groups that were
opened. A sensitivity pass, with the chosen layouts held fixed, changes the
node count, the block size, the replication factor, and fails node 0.

The same rules apply to every method. The codec model, the queries, and the
cluster do not change when the layout changes. For each seed the study JSON
records the shortlist: label, width, pack, duplication, trio co-location,
predicted objective, and measured objective.

## 5. Simulator results

Median ratio of the objective to the best pedigree-blind tiling, then to the
exhaustive winner of the fit mix. A ratio below 1 is a win against that
reference. The fit column is the mixture in Section 2. The other columns add
one full build to that workload's own queries, so they are not on the same
scale as the fit column.

| Method | Fit / generic | Family / generic | Extended / generic | Held-out family / generic | Scan / generic | Fit / oracle |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PedPack | 0.286 | 0.527 | 0.636 | 0.597 | 0.392 | 1.000 |
| No-pedigree ranker | 0.335 | 0.703 | 0.680 | 0.720 | 0.349 | 1.168 |
| Variant oracle | 0.335 | 0.703 | 0.680 | 0.720 | 0.349 | 1.168 |
| Family oracle | 0.286 | 0.527 | 0.636 | 0.597 | 0.392 | 1.000 |
| Fixed family heuristic | 0.286 | 0.530 | 0.637 | 0.597 | 0.392 | 1.007 |
| Generic oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 3.498 |

PedPack's fit choice, on every seed, is a nuclear-family layout with
group-major order and 16 tiles per file. Tile width is 400, 1,600, and 800
variants on seeds 7, 11, and 19. Trio co-location is 1. Duplication is 1.26,
1.21, and 1.27. The fit objectives are 3.224, 2.592, and 3.345. Each matches
the exhaustive winner and the family oracle.

The no-pedigree ranker selects, on every seed, the all-samples layout at
width 100 and 16 tiles per file. That is also the variant oracle. Its fit
objectives are 3.797, 3.026, and 3.705, which is 1.18×, 1.17×, and 1.11×
PedPack. Its shortlist is 24 variant layouts and 12 index-chunk layouts. It
never measures a family layout. PedPack's shortlist is 12 family, 20 variant,
and 4 pedigree layouts, and measurement then keeps the family layout.

On family-centric queries the byte gap is the point of the grouping. PedPack
reads 9.0×10⁵, 7.7×10⁵, and 1.1×10⁶ bytes. The variant layout reads 6.6×10⁶,
4.8×10⁶, and 6.0×10⁶, which is 7.3×, 6.2×, and 5.5× more. Held-out families,
drawn after the choice, show the same direction: 6.9×10⁵ versus 5.9×10⁶,
6.9×10⁵ versus 5.3×10⁶, and 6.1×10⁵ versus 4.7×10⁶ (8.6×, 7.7×, and 7.7×).

Scans still prefer one group. On seed 7 the family layout reads 2.0×10⁷
bytes in 0.191 s and the variant layout reads 1.7×10⁷ bytes in 0.162 s. The
ratio tracks the duplication. The 35% scan share is not large enough to
overturn the family queries once a family query no longer reads every sample.

The best index-chunk layout is 3.5× more expensive than PedPack on the fit
objective. Trio co-location is 0.15–0.21. On seed 7 it shuffles 1.5×10⁷ bytes
on the fit queries. The family layout and the variant layout shuffle none,
because every trio they are asked about sits inside a group they open.

The fixed heuristic is the same layout as PedPack on seed 7 (width 400). On
seeds 11 and 19 the search prefers a wider tile and the heuristic is 2.624
and 3.369 against PedPack's 2.592 and 3.345. The grouping, not the width
search, is the gap against the variant layout.

Agglomeration and the three-generation expansion do not beat nuclear families
on this mix: the pedigree oracle and the family oracle have the same ratios.
The encoder is a working kinship embedding — parent and child land closer
than unrelated founders in the unit test — and it does not change the layout
the cost model selects.

The cost network, scored on layouts outside its pilot, has seconds MAPE 0.89
and objective Spearman 0.91. The pedigree-aware formula on the same layouts
has MAPE 0.22 and Spearman 0.99. The network still lands near the family
oracle (fit / oracle 1.000 at the reported precision). It is not required
for the gap between PedPack and the blind ranker.

At replication 3, failing node 0 leaves the PedPack layouts with no
unavailable blocks. File contents were not repacked in that pass. The
sensitivity table is in `results/pedigree/study.json`.

## 6. What this does and does not show

The result that survives is the following. When queries name 5–15% of
families on a substantial variant window, and records are large relative to
the block, a ranker that knows which samples share a group selects nuclear
families, and a ranker that sees only group sizes selects an all-samples
layout. The family layout is the exhaustive winner of the 35/40/25 mix. It
reads several times fewer bytes on family queries, including families it was
not fit on, and it reads more on scans, in line with storing overlapping
generations twice. Pedigree-blind chunks are worse than both, because they
co-locate arbitrary samples and then shuffle the trios they split.

The same comparison, re-fit from this config, is in `results/stress/summary.md`
and `results/stress/boundary.png`. The reliable cases are frozen as
`configs/positive/select_05.json`, `select_15.json`, and `mix_50.json`.
`docs/boundary.md` states the generator settings on each side of the boundary.
Touching 5%, 15%, or 40% of families, and a fit mix that is half scans and
half family queries, still separates the two rankers on every seed. At 5%
the ablation's fit objective is 1.12–1.19× PedPack and it reads 12–16× more
bytes on held-out families. At 40% the objective ratio falls to 1.00–1.13
and the held-out byte ratio to 2.3–2.8×. The layouts are still different.

Three generations and a cross rate of 0.6 remove the gap. Duplication of the
nuclear layout rises to 1.33–1.43, the family oracle is more expensive than
the variant oracle, and both rankers select all-samples. A 20% scan share
restores a family choice on one of the three seeds only. A larger draw
(102–147 samples, 3,200 variants) is the same story: the seed with
duplication 1.45 selects all-samples; the other two still select families,
with fit ratios 1.018 and 1.001. The search is tracking that cost. It is not
a failure to notice the pedigree.

The same split survives a larger cohort. `configs/scale.json` draws 180
founder couples and 8,000 variants. Seeds 7, 11, and 19 have 1,175–1,207
people and nuclear duplication 1.28–1.30. PedPack again selects families and
matches the family oracle. The blind ranker selects all-samples. Its fit
objective is 1.20–1.28× PedPack, and it reads 16.3–17.6× more bytes on
held-out families. The family layout pays for that with build time (76–79 s
against 70–71 s) and storage (4.7–4.9 s against 4.2–4.3 s). Query time is
where it wins (4.4–5.9 s against 6.5–7.8 s). Scan bytes are 1.15–1.17×
higher for the family layout.

`configs/scale_3k.json` is the same regime with 500 founder couples. Seeds 7,
11, and 19 have 3,370, 3,361, and 3,328 people, nuclear duplication
1.293–1.298. PedPack selects `count:32` nuclear families and matches the
family oracle. The blind ranker selects all-samples. Its fit objective is
1.244×, 1.211×, and 1.208× PedPack. Held-out family bytes are 16.8×, 17.3×,
and 16.7×. Build and storage are higher for the family layout (about 218–220 s
and 13.4–13.6 s, against 197–199 s and 11.8–11.9 s). Query time is lower
(14.7–17.5 s against 20.0–23.5 s).

`configs/scale_overlap.json` is the closed side at the same variant count:
667–685 people, 3 generations, cross rate 0.6, family-oracle duplication
1.42–1.44. Both rankers select all-samples. The family oracle is more
expensive (seed 7: 218.7 against 173.1). Raising the sample count on that
3-generation draw did not reopen the gap. Duplication did.

The 1000 Genomes phase 3 pedigree is the same selector on a public graph.
`experiments/prepare_1kg_ped.py` keeps the 3,375 people who have both parents
or neither (636 children, 2,739 founders) from `20130606_g1k.ped`. Dosages
are still simulated. Nuclear duplication is 1.001. Seeds 7 and 11 have
PedPack on nuclear families at width 4,000 and the blind ranker on a
full-width family or agglomerated layout; the fit objectives differ by
0.15% and 0.007%. Seed 19 chooses the same agglomerated layout in both
rankers. All-samples remains the expensive layout: its fit objective is
1.42×, 1.43×, and 1.49× PedPack, and it reads 20.1×, 18.6×, and 19.9× more
bytes on held-out families. With essentially no overlapping membership, group
size is enough to reject all-samples. Pedigree identity is not doing extra
work. The run is `results/public_1kg/`.

An earlier configuration, 48 bytes per genotype and 256 KB blocks with a
25% scan share, did not produce this gap. Both rankers kept an all-samples
layout, and the numbers are in `results/small/`. That run is a different
workload. It is not a contradiction of the table above, and it is not
evidence that pedigree features help when the query still behaves like a
scan.

The result that does not survive is a claim that the neural ranker, or the
kinship encoder, or a three-generation overlap, is why. On this simulator
the formula ranks the pool better than the network, and nuclear families
match the best pedigree grouping in the pool. The fixed heuristic is close,
because the winning group structure is nuclear families; search still moves
the width.

A single-host Hadoop 3.2.1 cluster with one NameNode and three DataNodes
reproduced the byte ranking for the selective regimes and the tie for the
three-generation pedigrees. Details and the limits of that run are in
`results/hdfs/report.md`. It is one machine, the client is not a DataNode,
and there is still no remote-read fraction.

Hail 0.2.139 was run in Spark local mode on the same 1,207 × 8,000 matrix,
seed 7. Two comparisons, both local, neither a cluster.

The first (`results/hail/`) put the family layout in one exploded Table.
That object is the wrong encoding: 56–58 MB against a 3 MB MatrixTable,
while nuclear duplication is only 1.30. Query times sat under a second and
did not reproduce the 16× byte gap.

The second (`results/hail_columnar/`) uses the same integer MatrixTable for
both sides. The all-samples table is repartitioned into 16 locus partitions
(18.7 s to build, 3.0 MB). The family side is one MatrixTable per nuclear
family, 361 tables (937 s to build, 33.9 MB). Duplication is 1.30, so most
of the extra Hail bytes are per-table overhead. A held-out family query
returns the same 250,518 entries from both and takes 0.35–1.7 s on the
MatrixTable against 7.1–8.1 s on the 18 family tables. The scan takes 0.55 s
against 169 s. That run is a naive object model: one table per nuclear
family, 361 objects, where the search packed those groups 32 to a file.
It is not evidence about the layout. It is evidence that Hail plans each
table separately, and that cost dominates a 3 MB matrix.

The pack the search actually chose is measured in Parquet
(`results/parquet/`). Seed 7, Zstandard, row groups of 800 variants, 32
nuclear groups per file, 12 files. Build is 1.69 s and 4.54 MB for the
all-samples file, 1.28 s and 5.99 MB for the packed files (1.32×, next to
duplication 1.30). Held-out compressed bytes are 84.3 KB against 85.2 KB.
The packed read is slower (0.098 s against 0.037 s) because it opens 8–10
files. The scan reads 1.30× the bytes, which is the duplication. Column
projection already returns the named samples from the single file, so there
is no 16× byte gap left for packing to win.

`results/parquet_unified/` puts those same duplicated columns into one
file. On its own held-out draw the 12-file pack takes 0.106 s and 185,808
compressed bytes. The one-file pack, still row groups of 800 variants,
takes 0.063 s and the same 185,808 bytes. All-samples on that draw is
0.058 s and 181,150 bytes. One footer removes the file-open gap. The
selective bytes stay within about 3% of projecting the all-samples file,
which is the duplicated people the cover still reads, not a 16× tile.
A single row group over all 8,000 variants is faster still (0.031 s) and
reads more of each column (236,387 bytes), because the window can no
longer skip. On this matrix the row-group footer is a large part of a
4–6 MB file of one-byte alleles. The opaque-block gap is still the gap
of a reader that cannot name a column.

Nothing here is a multi-node HDFS run.
This host has 16 GB of RAM; four isolated virtual machines were not started.
The AWS, GCP, and Azure command-line tools are not installed on it, so an
EMR or Dataproc cluster was not started either. The block model has no
NameNode heap, no disk, and no scheduler noise. The
scaled cohorts run from about 1,200 people to about 3,370, at 8,000 variants, with a byte scale
factor so that a family tile spans several 64 KB blocks. They are not a
terabyte of 1000 Genomes. A paper that stops at this simulator, this
single-host Hadoop, and this local Hail timing should not claim a
distributed-systems speedup. The protocol those runs have to follow is in
`docs/evaluation.md`.

## 7. Related work

CoHadoop (Eltabakh et al., PVLDB 2011) co-locates HDFS blocks that share a
locator so a join can avoid a shuffle. PedPack also co-locates related
records. The unit is a genotype tile, the relatedness is the pedigree, and
the cost is allowed to refuse the co-location when scans and duplication
make it expensive. CoHadoop does not choose tile shape.

Trojan layouts (Jindal, Quiane-Ruiz, Dittrich, PVLDB 2011) pick a column
layout per HDFS block for a declared workload. The search here is that idea
on the sample axis of a matrix, with overlapping groups billed as extra
bytes rather than as a join key.

A Hail MatrixTable is partitioned by locus interval. Filtering to a few
samples still reads every sample in the opened partitions. That is the
all-samples layout in this paper. ADAM (Massie et al., UCB/EECS-2013-207)
and Glow put reads and variants on Spark in the same locus-major shape.
They standardize the schema. They do not search a pedigree cost.

TileDB (Papadopoulos et al., PVLDB 2016) stores a dense array as tiles whose
shape is chosen by the user. PedPack is a policy for the sample side of a
dosage matrix when some queries name families and the groups overlap. The
array store does not supply that policy.

## 8. Where the code matches the text

`src/pedpack/study.py` is the protocol. `cost.py` is the formula, including
the blind ranker. `boundary.py` is the two-layout tradeoff.
`simulator.py` is the measurement. `workloads.py` builds the family-centric,
extended, and held-out queries. `selector.py` is the budget. The encoder and
the cost network in `models.py` are not the selector.
`configs/pedigree_gap.json` is the locked small-cohort run. `configs/scale.json`
is the 1,200-person run. `configs/scale_3k.json` is the 3,300-person run.
`configs/public_1kg.json` is the 1000 Genomes pedigree with simulated dosages.
`results/pedigree/report.md` is the table this draft was written from.
`results/compression/report.md` is the Zstandard check of the kinship term.
