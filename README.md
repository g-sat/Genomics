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
debug the cost model. It is not a stand-in for a multi-node HDFS run or for
Hail. `docs/evaluation.md` says which of those runs are still required before
the systems claim can be made at full strength.

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
python experiments/run_study.py --config configs/small.json --out results/small
pytest
```

`configs/tiny.json` is the smoke study used by the tests. `configs/small.json`
is the simulator comparison behind `docs/manuscript.md`. The current run is
summarized in `results/small/report.md`.

On that run, with a 32-query horizon and a 25% scan share, the selector
matches an exhaustive search and picks an all-samples layout. That layout's
mixed objective is about 0.41× the best pedigree-blind tiling. A fixed
nuclear-family packing is worse on the same mix because multi-generation
families duplicate people. Retarget the selector at family-scoped queries and
it moves to family groups. The neural ranker is less accurate than the
closed-form model on this simulator; the gain is the search space and the
objective, not the network. Details and the limits of the simulator are in
the manuscript.
