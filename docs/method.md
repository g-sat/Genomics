# Method

PedPack lays out a dosage matrix of `M` variants by `N` samples. The pedigree
is part of the physical design, not a join done later.

## Layout space

A logical tile is one sample group over a contiguous variant span. Three
decisions define a layout:

- the variant width
- the sample groups
- the packing of tiles into files (`naive` is one tile per file, `count:k`
  concatenates `k` tiles, `bytes:B` fills files up to `B` bytes)

Write order is either group-major (a family's tiles are adjacent) or
locus-major (a variant span's groups are adjacent). The planner reads whole
tiles. A sample that appears in several groups is stored more than once;
storage pays for every copy, and the query planner reads one covering set of
groups.

Groups in the search pool:

| Family tag | How the groups are built |
| --- | --- |
| `variant` | One group containing every sample |
| `generic` | Index-contiguous chunks. Columns are scrambled so this is not secretly a family layout |
| `family` | Nuclear family: two parents and their children |
| `pedigree` | Size-capped kinship agglomeration, and optional three-generation groups that overlap |
| `heuristic` | Nuclear families at one fixed width and one fixed pack. No search |

Agglomeration starts from individuals and merges the pair with the highest
score that does not exceed a size cap. The score is mean cross-kinship plus
`λ` times the cosine of the two clusters' mean encoder embeddings. `λ = 0` is
kinship only.

## Objective

For a layout `L` and a fit mixture with shares `α_scan`, `α_family`, and
`α_extended`, amortized over a horizon of `H` queries:

```
C(L) = T_build(L) + ω · S(L)
     + H · (α_scan · mean T_scan
            + α_family · mean T_family
            + α_extended · mean T_extended)
```

The pedigree-gap configuration uses `α_scan = 0.35`, `α_family = 0.40`,
`α_extended = 0.25`, and `H = 32`. A raw sum of the query list is not used for
selection: one scan reads much more than a family-scoped query, so an
unweighted sum is a scan workload. Shares of empty buckets are renormalized
onto the buckets that have queries.

`S(L)` is the replicated size. `ω` prices that size as `storage_scan_equivalents`
reference full scans of the all-samples layout. The default is `0.25`. Build,
storage, and query time are seconds.

Build time reads the source matrix once, shuffles a fraction `(1 − trio
co-location)` of it, writes the primary copy locally, and writes extra
replicas at remote bandwidth.

Tile bytes are

```
raw = samples × variants × bytes_per_genotype
factor = clip(compress_base × (1 − kinship_alpha × mean_within_kinship), 0.08, 1)
```

`factor` stands in for a codec that gets smaller on related haplotypes. It is
not a PBWT implementation. Replacing `payload_bytes` is the intended seam.

## Trio co-location and packing

A trio is co-located when one sample group contains the father, the mother,
and the child. The layout rate is that fraction over every trio in the
pedigree. The query rate is that fraction over the trios named in the query.
Variant span is not part of the test: every group is tiled across the variant
axis. An all-samples group has co-location 1 because every trio is inside it.

Packing is independent of grouping:

- `naive` writes one tile per file
- `count:k` concatenates `k` tiles in write order
- `bytes:B` starts a new file when the next tile would pass `B` bytes

Group-major order keeps one group's tiles together, so a family query tends
to open fewer files. Locus-major order keeps one variant span's groups
together, which helps a scan and hurts a family query.

## What the formula is allowed to see

`cost.analytical_query` ranks a layout without walking replicas. With
`pedigree_aware=True` it opens the true covering groups, uses the
kinship-compressed tile sizes, rounds each touched file up to a whole block
(capped at the file size), and charges shuffle from this query's trios:

```
shuffle = (1 − query_trio_colocation) × n_trios × span × 3 × bytes_per_genotype
```

The remote fraction is still not task placement. It is

```
min(0.85, 0.05 × max(0, 1 − replication / n_nodes) + 0.25 × (files − 1) / files)
```

so a one-file layout stays near the replication baseline and a many-file
layout pays more remote traffic.

With `pedigree_aware=False` the ranker keeps group sizes and packing and
drops identity. A group of size `s` is opened with probability
`1 − (1 − q/N)^s`, tiles are sized with kinship 0, and co-location is
`mean group size / N`. Two layouts with the same size histogram get the same
selective-byte estimate. That is the ranker used by
`pedpack_no_pedigree_features`.

## Simulator

`simulator.measure` is the measurement the formula is scored against.

- One task is opened per file that the query touches.
- The task runs on the live node holding the most replicas of the touched
  blocks. Other blocks of that file are remote.
- Replica `r` of block `b` in file `f` is node `(h + r) mod n_nodes`, with
  `h = (f · 10007 + b · 97) mod n_nodes`.
- A block whose replicas are all on failed nodes is unavailable and is charged
  an extra remote read.
- A trio whose three members do not share a selected group contributes
  `3 × variant_span × bytes_per_genotype` of shuffle.

The pedigree-aware formula and the simulator use the same greedy group cover.
They can still disagree on remote bytes, because the formula does not place
the task. The blind formula does not use that cover.

## Kinship encoder

Each person is a node. Edges are parent–child and mate links, stored in both
directions. Node features are founder status, sex, log child count, parent
count, and mean kinship. Two rounds of mean aggregation produce an embedding.
Training minimizes the squared error between cosine similarity and kinship
divided by 1/2, clipped to `[0, 1]`. Query times are not in this loss. The
encoder runs before any layout is measured.

## Cost network

The network is a two-hidden-layer perceptron. Inputs are the columns in
`features.FEATURE_NAMES`: tile width, group-size statistics, trio
co-location, parent–child cut rate, within-group kinship, duplication, file
and block counts, pack and order flags, replication, and query shape
(sample fraction, variant fraction, query-specific trio co-location, read
amplification). Four outputs are log bytes, remote fraction, log shuffle, and
log milliseconds.

Rows and targets are standardized on the training split. Early stopping uses
a held-out slice of those rows. At prediction time the outputs are mapped
back to seconds. Build and storage stay in the closed form, because both are
known once the files are described.

Pedigree columns can be zeroed. That ablation keeps the same candidate pool
and the same pilot, and removes relatedness from the ranker. A second
ablation keeps the ranker and deletes every pedigree grouping from the pool,
leaving variant and index-chunk layouts.

## Search under a measurement budget

Let `B` be the number of layouts the selector may measure on the fit mix.
`pedpack` ranks every layout with the pedigree-aware formula and measures the
top `B`. The returned layout is the best among those measurements, not the
formula's first place if measurement disagrees.

`pedpack_no_pedigree_features` is the same budget with the blind formula.
`analytical_budget` is the same search as `pedpack`. `random_budget` measures
`B` random layouts.

The reported studies do not train a cost network. The code can still fit one
(`pedpack_neural`) on a stratified pilot, but that network is not the
selector and it is not part of the pedigree claim. A report that skips it
says so, instead of printing a missing error as a score.

`family_versus_variant` in `boundary.py` is the same objective reduced to
all-samples versus nuclear families. It is the map of duplication against
the fraction of families a query names.

Oracles measure every layout in a family. They are upper bounds.
`full_oracle` is the best layout in the whole pool on the fit mix.
`pedigree_oracle` is restricted to family, agglomerative, and heuristic
layouts. `variant_oracle` is restricted to all-samples layouts.
`generic_oracle` is restricted to index chunks. `heuristic_fixed` is one
nuclear layout at a fixed width and pack.

The fit mix is the only workload used to choose `pedpack`. Scan,
family-selective, extended, and held-out lists are scored after the choice is
frozen. `analysis.best_on_heldout_family` is allowed to see the held-out
family queries. It is an upper bound, not a method.

Every measured shortlist row is stored under `decisions` in the study JSON:
label, family, width, pack, order, file count, duplication, trio
co-location, predicted objective, and measured objective.
