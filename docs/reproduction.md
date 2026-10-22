# Reproduction

Requires Python 3.11 or newer, NumPy, and PyTorch. Matplotlib is only used to
write `fit_ratios.png`.

```
pip install -e ".[dev,plot]"
pytest
python experiments/run_study.py --config configs/pedigree_gap.json --out results/pedigree
python experiments/run_stress.py
```

That run writes `results/pedigree/study.json`, `results/pedigree/report.md`,
and, if Matplotlib is installed, `results/pedigree/fit_ratios.png`.
`configs/small.json` is the earlier cluster scale (48 bytes per genotype,
256 KiB blocks). `results/small/` was produced before the family-selective
workloads and should not be compared number-for-number with `results/pedigree/`.

`configs/tiny.json` is the configuration the tests execute. It is there to
check the protocol (budget, coverage, Mendelian consistency), not to estimate
a workload.

## Using a real pedigree

PED files are PLINK-style: family, individual, father, mother, sex,
phenotype. Phenotype is ignored. Parent id `0` means a founder. A parent who
is absent from the file is treated as unknown.

```python
from pedpack.genotype import cohort_from_ped
import numpy as np

cohort = cohort_from_ped("cohort.ped", n_variants=1600, rng=np.random.default_rng(0))
```

That call keeps the pedigree and simulates dosages. To use real hard calls,
load a `(variants, samples)` array of values in `{0, 1, 2}` aligned to the
PED order before column permutation:

```python
from pedpack.genotype import load_dosage_npy
from pedpack.pedigree import parse_ped

pedigree = parse_ped("cohort.ped")
dosages = load_dosage_npy("dosages.npy", pedigree)
```

Column order in the study is then scrambled so that a contiguous chunk is not
accidentally a family. Apply the same permutation to the dosage matrix that
`build_cohort` applies to simulated dosages.

## Larger cohorts

```bash
python experiments/run_study.py --config configs/scale.json --out results/scale
python experiments/run_study.py --config configs/scale_overlap.json --out results/scale_overlap
```

The first is 180 founder couples, 2 generations, 8,000 variants, family
queries at 5% of families. The second is 80 founder couples, 3 generations,
cross rate 0.6. Both use seeds 7, 11, and 19. Agglomeration at this size is
the average-linkage update in `grouping.py`; it is the same merge rule as
the small pedigrees.

Hail, if installed, times a variant-keyed MatrixTable against a family-keyed
table on the first config:

```bash
python experiments/hail_baseline.py --config configs/scale.json --out results/hail
```

The run in this repository used Hail 0.2.139 in Spark local mode inside WSL.
That is one virtual machine, not a Hail cluster. One MatrixTable per nuclear
family is the naive object counted in `results/hail_columnar/`. The pack of
32 groups per file, as Parquet, is:

```bash
python experiments/measure_parquet.py
python experiments/measure_parquet_unified.py
```

The larger simulator cohort and the Zstandard kinship check are:

```bash
python experiments/run_study.py --config configs/scale_3k.json --out results/scale_3k
python experiments/measure_compression.py
python experiments/prepare_1kg_ped.py
python experiments/run_study.py --config configs/public_1kg.json --out results/public_1kg
```

## Seeds

`configs/pedigree_gap.json` and `configs/small.json` use seeds 7, 11, and 19. The simulator is
deterministic given the seed. A second run of the same config produces the
same objectives. Variance in the report is variance across pedigrees and
query draws.

## What not to change between methods

The study already holds the cohort, the query definitions, the byte model,
the execution rules, and the cluster fixed while it swaps layouts. When you
move a winner onto HDFS, keep the codec, the replication factor, and the
query answers fixed the same way. Compare bytes read and shuffle from the
storage system, and check answers against the full-matrix Mendelian counts.
