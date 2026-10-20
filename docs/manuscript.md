# Pedigree-Aware Physical Design for Distributed Genomic Storage

## Abstract

We choose tile width, sample grouping, and file packing for a genotype matrix
whose queries mix full scans with family-scoped hereditary lookups. The
objective adds build time, replicated storage, and expected query time. Sample
groups may ignore the pedigree, follow nuclear families, or follow a kinship
encoder. A closed-form cost model and a small neural network both rank
candidate layouts; the layout we keep is the one a block simulator measures as
cheapest, inside a fixed measurement budget.

On three synthetic cohorts (58–76 samples, 1,600 variants, byte costs scaled
to a 256 KB block) a 32-query horizon with a 25% scan share selects an
all-samples layout. That layout matches an exhaustive search of 240 candidates.
Its mixed objective is 0.41× the best pedigree-blind two-dimensional tiling.
The blind tiling splits trios and pays application shuffle. A fixed nuclear-family
heuristic is also worse on this mix: people who are both children and parents
are stored in two families (duplication 1.21–1.27), and scans read the copies.
Retargeting the same selector at family-scoped queries switches it to family
groups and, on two of the three seeds, cuts selective bytes from about
7×10⁵ to about 1–2×10⁵. A message-passing kinship encoder did not beat nuclear
families. The neural ranker was less accurate than the closed-form model
(seconds MAPE 0.38 versus 0.24; objective Spearman 0.71 versus 0.98) and was
not where the gain came from. These numbers are from a block simulator. They
are not HDFS or Hail measurements.

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

The reason to search this space, rather than to declare that related samples
belong together, is that "together" has two cheap implementations. One tile
of everyone is together. One tile per nuclear family is together for that
family and expensive for a scan. Which one wins depends on the fraction of
scans, the duplication, and the file-open cost. The selector exists to make
that call from a cost, not from a slogan.

## 2. Cost

For a layout L and a horizon of H queries,

```
C(L) = T_build(L) + ω S(L) + H (α mean T_scan + (1 − α) mean T_selective)
```

S(L) is the replicated size. ω is set so that storing the reference
all-samples layout costs 0.25 reference scans. α is 0.25 and H is 32 in the
reported run. Build, storage, and query time are all in seconds.

An unweighted sum of a query list is the wrong objective here. A scan reads
the matrix; a family-scoped query reads a few percent of it. Four scans plus
ten selective queries is a scan benchmark. The mixture makes the scan share
a parameter instead of an accident of list length.

Tile bytes are the raw record size times a compression factor

```
clip(0.45 × (1 − 0.55 × mean within-group kinship), 0.08, 1)
```

That factor is a stand-in for a codec that shrinks on related haplotypes. It
is not a measured PBWT. Build time reads the source once, shuffles a fraction
(1 − trio co-location) of it, writes the primary copy, and writes extra
replicas at remote bandwidth.

The closed-form query model is missing two effects the simulator includes.
Bytes are the overlapping fraction of a tile, with no round-up onto the block
that is actually read. The remote fraction is a function of replication and
node count, not of which node the task was placed on. Shuffle uses the
layout-wide trio co-location rate, so a query against one intact family is
over-charged on a layout that splits other families. Those gaps are the
reason a learned ranker could, in principle, beat the formula.

## 3. What is learned

The kinship encoder is a two-layer mean aggregator on the pedigree. Edges are
parent–child and mate links. It is trained so that cosine similarity matches
kinship / ½, clipped to [0, 1]. Query times are not in that loss. Agglomeration
may add λ times the embedding cosine to the kinship score. λ = 0 is in the
pool as well.

The cost network is a two-hidden-layer perceptron. It maps layout and query
features to log bytes, remote fraction, log shuffle, and log milliseconds.
Features include tile width, group size, trio co-location, the parent–child
cut rate, within-group kinship, duplication, file and block counts, and the
query's own sample and variant fractions. They do not include simulator
outputs. Targets are standardized. Early stopping uses a held-out slice of
the training rows.

The selector may measure B = 36 layouts out of 240. It measures a stratified
pilot of about half of that budget, fits the network, and spends the rest on
the layouts the network ranks highest. The returned layout is the cheapest
among the ones that were measured. An analytical budget spends all 36
measurements on the formula's top list. A random budget measures 36 random
layouts. Oracles measure every layout in a family. They are upper bounds.
`pedpack_no_pedigree_features` zeros the pedigree columns and keeps the same
pilot. `pedpack_generic_pool` deletes every pedigree grouping from the pool.

## 4. Protocol

Cohorts are simulated: 10 founder couples, two generations of children, a
25% chance that a couple spans birth families. Dosages are transmitted
without Mendelian error; the study aborts if that check fails. Columns are
permuted so that consecutive samples are not a family. Each seed draws its
own pedigree and its own queries. Seeds are 7, 11, and 19.

Fit queries are large ranges over every sample, plus family-scoped ranges on
training families. Selective, scan, held-out, and short-window workloads are
scored after the layout is chosen. Held-out queries use families that were
withheld from the fit list. The window queries are short spans on those
families. They imitate the access pattern of a compound-heterozygote search.
They do not call variants.

The simulator places replicas on 8 nodes at replication 3, reads a block from
the node that holds the most replicas of the file's touched blocks, and
charges shuffle when a queried trio is split across the groups that were
opened. A sensitivity pass, with the chosen layouts held fixed, changes the
node count, the block size, the replication factor, and fails node 0.

The same rules apply to every method. The codec model, the queries, and the
cluster do not change when the layout changes.

## 5. Simulator results

Median ratio of the objective to the best pedigree-blind tiling, then to the
exhaustive winner of the fit mix. A ratio below 1 is a win against that
reference. The fit column is the mixture in Section 2. The other columns add
one full build to that workload's own queries, so they are not on the same
scale as the fit column.

| Method | Fit / generic | Scan / generic | Selective / generic | Held-out / generic | Fit / oracle |
| --- | ---: | ---: | ---: | ---: | ---: |
| PedPack (fit mix) | 0.409 | 0.405 | 0.752 | 0.728 | 1.000 |
| Analytical budget | 0.409 | 0.405 | 0.752 | 0.728 | 1.000 |
| Variant oracle | 0.409 | 0.405 | 0.752 | 0.728 | 1.000 |
| Family oracle | 0.464 | 0.521 | 0.680 | 0.594 | 1.004 |
| Fixed family heuristic | 0.491 | 0.561 | 0.677 | 0.588 | 1.069 |
| Pedigree oracle | 0.464 | 0.521 | 0.680 | 0.594 | 1.004 |
| Generic oracle | 1.000 | 1.000 | 1.000 | 1.000 | 2.445 |
| PedPack tuned on selective queries | 0.515 | 0.563 | 0.672 | 0.578 | 1.310 |

PedPack's fit choice, on every seed, is an all-samples layout: one group,
duplication 1, trio co-location 1, one or two files, tile width 100 or 800
variants. It matches the exhaustive winner. The best index-chunk layout is
2.4× more expensive on the same objective. Its trio co-location is 0.15–0.21.
On seed 7 that layout shuffles 1.9×10⁶ bytes on the fit queries; the
all-samples layout shuffles none.

The fixed family heuristic loses to that all-samples layout on the mix and
on scans. Duplication is 1.21, 1.26, and 1.27 on the three seeds. A scan of
the family layout on seed 7 reads 2.9×10⁶ bytes against 1.7×10⁶ for the
all-samples layout. The extra bytes are the overlapping generation, not a
failure of co-location: trio co-location on the family layout is also 1.

Family groups win when the objective ignores scans. On seed 7 the
selective-tuned selector picks a family layout and drops selective bytes from
6.4×10⁵ to 2.2×10⁵, and held-out bytes from 1.1×10⁶ to 3.1×10⁵. Held-out
families were not in the selective training queries. On seed 11 the
selective-tuned layout is again a family layout (naive packing, 16 files):
selective seconds fall from 0.0085 to 0.0046, and scan seconds rise from
0.016 to 0.077. On seed 19 the selective-tuned search returns a kinship
grouping with trio co-location 0.79, which is worse on selective seconds
(0.0072) than the nuclear-family oracle (0.0054). That is selector regret,
and it lines up with the ranker's accuracy.

The cost network, scored on layouts it did not train on, has seconds MAPE
0.38 and objective Spearman 0.71. The formula on the same layouts has MAPE
0.24 and Spearman 0.98. Hiding pedigree features does not change the fit-mix
choice, because the winning layout is an all-samples layout and the formula
already ranks that plateau correctly. Random measurement of 36 layouts also
hits the plateau. The budget is easy when many variant layouts are nearly
tied. It is not easy on the selective retune, where the network missed the
best family layout on one seed.

Agglomeration and the three-generation expansion do not beat nuclear families
on the fit mix: the pedigree oracle and the family oracle have the same
median ratios. The encoder is a working kinship embedding — parent and child
land closer than unrelated founders in the unit test — and it does not change
the layout the cost model selects.

At replication 3, failing node 0 does not change bytes or the remote fraction
of the fit-mix layout on seed 7. At replication 1 the remote fraction on that
layout rises from 0.12 to 0.62. A 64 KB block raises the remote fraction to
0.55. A 512 KB block drops it to 0 and raises bytes read, because a narrow
query pulls a larger block. File contents were not repacked for that pass.
The layout was chosen at 256 KB and replication 3.

## 6. What this does and does not show

The result that survives is the following. For a mix that still contains a
quarter of its queries as scans, the cheapest layout co-locates every sample
in a few large files. Family packing is the wrong default on that mix,
because multi-generation pedigrees duplicate people. Family packing is the
right default when the queries are family-scoped, including on families the
selector did not fit, and the same search procedure will move there if that
is the objective. Pedigree-blind chunks are worse than both, because they
co-locate arbitrary samples and then shuffle the trios they split.

The result that does not survive is a claim that the neural ranker, or the
kinship encoder, is why. On this simulator the formula ranks the pool better.
The encoder does not find a grouping that nuclear families miss.

Nothing here is a multi-node HDFS run, a Hail comparison, or a measurement of
a real codec. The block model has no NameNode heap, no disk, and no scheduler
noise. The cohorts are tens of samples, with a byte scale factor so that
tiles cross a 256 KB block; they are not a terabyte of 1000 Genomes. A paper
that stops at this simulator should not claim a distributed-systems speedup.
The protocol those runs have to follow is in `docs/evaluation.md`: the same
queries, the generic two-dimensional baseline, the fixed family heuristic,
build and storage in the objective, and a held-out family set. Hail, if it
is compared, has to be run. This repository does not contain a Hail number.

## 7. Where the code matches the text

`src/pedpack/study.py` is the protocol. `cost.py` is the formula. `simulator.py`
is the measurement. `models.py` is the encoder and the ranker. `selector.py`
is the budget. `configs/small.json` is the run behind Section 5.
`results/small/report.md` is the table this draft was written from.
