# Reproduction

Requires Python 3.11 or newer, NumPy, and PyTorch. Matplotlib is only used to
write `fit_ratios.png`.

```
pip install -e ".[dev,plot]"
pytest
python experiments/run_study.py --config configs/pedigree_gap.json --out results/pedigree
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
