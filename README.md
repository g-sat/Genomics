# PedPack

PedPack chooses how a genotype matrix and its pedigree are cut into tiles and
packed into block-storage files. The objective is one number, in seconds:

```
build + (replicated bytes × storage price) + expected query time
```

Sample groups can follow nuclear families, kinship agglomeration, or an
overlapping three-generation expansion. A two-layer message-passing encoder,
trained only to reconstruct kinship, can bias the agglomeration. A small cost
network, trained on a pilot of measured layouts, ranks the rest. The selector
remeasures a shortlist and keeps the best measured layout. It does not trust
the network for the final choice.

The measurement in this repository is a block simulator: file offsets, block
boundaries, replica placement, and a failed DataNode. It is the right place to
debug the cost model. It is not a stand-in for a multi-node HDFS run.
`results/hail/` is Hail 0.2.139 in local mode on the scaled matrices, and
`docs/evaluation.md` says which runs are still required before the systems
claim can be made at full strength.

## Layout

```
configs/          study settings
docs/             method, protocol, reproduction
experiments/      study driver
src/pedpack/      library
tests/            invariants
results/          generated reports (not source)
```

| Module | Role |
| --- | --- |
| `pedigree.py`, `genotype.py` | PED I/O, kinship, Mendelian dosages |
| `grouping.py` | Nuclear, generic, agglomerative, overlapping groups |
| `layout.py`, `manifest.py` | Tiles, packing, writer manifest |
| `cost.py` | Closed-form cost, deliberately short of the simulator |
| `simulator.py` | Block placement and query measurement |
| `models.py` | Kinship encoder and cost network |
| `selector.py` | Budgeted search |
| `study.py` | Baselines, ablations, held-out queries |

## Run

```
pip install -e ".[dev,plot]"
python experiments/run_study.py --config configs/pedigree_gap.json --out results/pedigree
pytest
```

`configs/tiny.json` is the smoke study used by the tests.
`configs/pedigree_gap.json` is the simulator comparison behind
`docs/manuscript.md`. The current run is summarized in
`results/pedigree/report.md`. `results/small/` is an earlier workload
(48 bytes per genotype, 256 KB blocks, 25% scan share) on which an
all-samples layout won. Do not mix the two tables.

On the pedigree-gap run, with a 32-query horizon split 35% scans, 40%
family-centric queries, and 25% extended-pedigree queries, PedPack matches
an exhaustive search and picks nuclear families on every seed. The same
budget with a ranker that cannot see who shares a group picks the all-samples
layout, and that layout reads about 5–9× more bytes on family and held-out
family queries. Scans still cost more on the family layout, because people
in two generations are stored twice. The neural ranker is less accurate than
the closed-form model on this simulator. Details and the limits of the
simulator are in the manuscript.
