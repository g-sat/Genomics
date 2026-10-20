# Pedigree-Aware Physical Design for Distributed Genomic Storage

## Abstract

We choose tile width, sample grouping, and file packing for a genotype matrix
whose queries mix full scans with family-scoped hereditary lookups. The
objective adds build time, replicated storage, and expected query time. The
selector that we call PedPack ranks candidate layouts with a closed-form cost
model that can see which samples share a group, then keeps the cheapest layout
a block simulator measures inside a fixed budget. A second ranker sees the
same group sizes and the same packing, and does not see who is related to
whom.

On three synthetic cohorts (58–76 samples, 1,600 variants, 400 bytes per
genotype, 64 KB blocks) a 32-query horizon with a 35% scan share selects a
nuclear-family layout on every seed. That layout matches an exhaustive search
of 240 candidates. The ranker that cannot see pedigree identity selects the
all-samples layout on every seed, the same layout the variant oracle selects.
The all-samples fit objective is 1.11–1.18× the family layout. On family-scoped
queries it reads 5.5–7.3× fewer bytes, and on held-out families, which were
not used to choose the layout, 7.6–8.6× fewer. Scans go the other way: the
family layout stores each multi-generation person twice (duplication
1.21–1.27) and reads about 1.2× more on a scan. Index-chunk layouts are worse
than both, because they split trios and pay shuffle. A neural ranker was less
accurate than the formula (seconds MAPE 0.89 versus 0.22; objective Spearman
0.91 versus 0.99) and is not the method behind the gap. These numbers are
from a block simulator. They are not HDFS or Hail measurements.

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
on the chosen layout). Rank correlation on layouts held out of the neural
pilot is 0.99, which is what the search uses.

## 3. What is learned, and what actually searches

The kinship encoder is a two-layer mean aggregator on the pedigree. Edges are
parent–child and mate links. It is trained so that cosine similarity matches
kinship / ½, clipped to [0, 1]. Query times are not in that loss. Agglomeration
may add λ times the embedding cosine to the kinship score. λ = 0 is in the
pool as well.

The cost network is a two-hidden-layer perceptron. It maps layout and query
features to log bytes, remote fraction, log shuffle, and log milliseconds.
It is a secondary experiment. PedPack itself does not use it.

PedPack ranks all 240 layouts with the pedigree-aware formula and measures
the top B = 36. The returned layout is the cheapest among those
measurements. `pedpack_no_pedigree_features` is the same budget with the
blind formula. `analytical_budget` is the same search as PedPack.
`random_budget` measures 36 random layouts. Oracles measure every layout in
a family. They are upper bounds. `pedpack_neural` and `pedpack_neural_blind`
are the network with and without pedigree feature columns. They are not the
pedigree claim.

## 4. Protocol

Cohorts are simulated: 10 founder couples, two generations of children, a
25% chance that a couple spans birth families. Dosages are transmitted
without Mendelian error; the study aborts if that check fails. Columns are
permuted so that consecutive samples are not a family. Each seed draws its
own pedigree and its own queries. Seeds are 7, 11, and 19.

A family-centric query names 5–15% of the nuclear families and a window
covering 20–40% of the variant axis. An extended query starts from one
nuclear family, adds the parents of its members, and with probability one
half adds a second family, so the request is not a clean independent trio.
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

Nothing here is a multi-node HDFS run, a Hail comparison, or a measurement of
a real codec. The block model has no NameNode heap, no disk, and no scheduler
noise. The cohorts are tens of samples, with a byte scale factor so that a
family tile spans several 64 KB blocks; they are not a terabyte of 1000
Genomes. A paper that stops at this simulator should not claim a
distributed-systems speedup. The protocol those runs have to follow is in
`docs/evaluation.md`: the same queries, the generic two-dimensional baseline,
the fixed family heuristic, the blind ranker, build and storage in the
objective, and a held-out family set. Hail, if it is compared, has to be
run. This repository does not contain a Hail number.

## 7. Where the code matches the text

`src/pedpack/study.py` is the protocol. `cost.py` is the formula, including
the blind ranker. `simulator.py` is the measurement. `workloads.py` builds
the family-centric, extended, and held-out queries. `models.py` is the
encoder and the secondary network. `selector.py` is the budget.
`configs/pedigree_gap.json` is the run behind Section 5.
`results/pedigree/report.md` is the table this draft was written from.
