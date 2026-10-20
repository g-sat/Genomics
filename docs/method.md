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

For a layout `L` and a fit mixture with scan share `α` and selective share `1 − α`,
amortized over a horizon of `H` queries:

```
C(L) = T_build(L) + ω · S(L) + H · (α · mean T_scan + (1 − α) · mean T_selective)
```

The default mixture is `α = 0.25`, `H = 32`. A raw sum of the query list is not
used for selection: one scan reads orders of magnitude more than a family-scoped
query, so an unweighted sum is a scan workload with a selective footnote.
Scan-only and selective-only workloads are still reported on their own, and
`pedpack_tuned_on_scan` re-runs the search with the scan queries as the objective.

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

## What the formula is allowed to see

`cost.analytical_query` ranks a layout without walking replicas. It is wrong
on purpose in two places that matter on block storage:

- bytes are the overlapping fraction of each tile, with no round-up to the
  block that is actually read
- the remote fraction is `0.05 × max(0, 1 − replication / n_nodes)`, not the
  task placement

Shuffle uses the layout-wide trio co-location rate, not the trios in the
query. A selective query against one intact family is therefore over-charged
on a layout that splits other families.

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

The same greedy group cover is used by the formula and the simulator, so they
disagree about bytes and placement, not about which groups are opened.

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
The neural selector measures a stratified pilot of about `B/2` layouts
(variant, generic, and pedigree buckets), fits the network, predicts the
other layouts, and measures the best remainder until `B` measurements are
used. The returned layout is the best among the ones it actually measured.

Two budgeted controls use the same `B`:

- analytical: measure the `B` layouts the formula ranks first
- random: measure `B` random layouts

Oracles measure every layout in a family. They are upper bounds, labeled as
such. `full_oracle` is the best layout in the whole pool on the fit mix.
`pedigree_oracle` is restricted to family, agglomerative, and heuristic
layouts. `generic_oracle` is restricted to index chunks. The gap between
those two oracles is the value of pedigree groupings when search is free.
The gap between `pedpack` and `full_oracle` is the regret of the budgeted
search. The gap between `pedpack` and `pedpack_generic_pool` is the value of
being allowed to propose pedigree groups under the same budget.

The fit mix is the only workload used to choose `pedpack`. Scan, selective,
held-out families, and short windows are scored after the choice is frozen.
`pedpack_tuned_on_scan` repeats the search with the scan queries as the
objective, which shows whether the operating point moves when the mix changes.
`analysis.best_on_heldout` is allowed to see the held-out queries. It is an
upper bound, not a method.
