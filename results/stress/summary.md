# Stress tests around the locked pedigree-gap regime

Each scenario re-runs the search. Layouts are chosen on that scenario's
fit mix, then scored on held-out families. Numbers are from the block
simulator. They are not HDFS measurements.

Base config: seeds [7, 11, 19], horizon 32, mixture {'scan': 0.35, 'family': 0.4, 'extended': 0.25}, block 65536, bytes per genotype 400.

| Scenario | Samples | PedPack | No-pedigree | Blind / PedPack | Same layout? | Held-out byte ratio | Gap |
| --- | ---: | --- | --- | ---: | --- | ---: | --- |
| select_05 | 76 | family | variant | 1.185 | False | 13.641 | yes |
| select_05 | 58 | family | variant | 1.193 | False | 12.029 | yes |
| select_05 | 67 | family | variant | 1.121 | False | 15.719 | yes |
| select_15 | 76 | family | variant | 1.141 | False | 5.780 | yes |
| select_15 | 58 | family | variant | 1.211 | False | 5.963 | yes |
| select_15 | 67 | family | variant | 1.107 | False | 5.374 | yes |
| select_40 | 76 | family | variant | 1.066 | False | 2.528 | yes |
| select_40 | 58 | family | variant | 1.130 | False | 2.315 | yes |
| select_40 | 67 | family | variant | 1.004 | False | 2.768 | yes |
| mix_50 | 76 | family | variant | 1.097 | False | 6.944 | yes |
| mix_50 | 58 | family | variant | 1.088 | False | 5.000 | yes |
| mix_50 | 67 | family | variant | 1.040 | False | 8.828 | yes |
| overlap | 104 | variant | variant | 1.000 | False | 1.000 | no |
| overlap | 73 | variant | variant | 1.000 | True | 1.000 | no |
| overlap | 97 | variant | variant | 1.000 | True | 1.000 | no |
| scale | 147 | variant | variant | 1.000 | True | 1.000 | no |
| scale | 102 | family | variant | 1.018 | False | 7.535 | yes |
| scale | 147 | family | variant | 1.001 | False | 6.827 | yes |

## What each scenario changes

- `select_05`: 5% of families. Family-selective queries touch 5.1% of samples on seed 7. Gap on every seed: yes.
- `select_15`: 15% of families. Family-selective queries touch 19.1% of samples on seed 7. Gap on every seed: yes.
- `select_40`: 40% of families. Family-selective queries touch 47.2% of samples on seed 7. Gap on every seed: yes.
- `mix_50`: fit mix 50% scan / 50% family. Family-selective queries touch 10.7% of samples on seed 7. Gap on every seed: yes.
- `overlap`: 3 generations, cross rate 0.6. Family-selective queries touch 17.7% of samples on seed 7. Gap on every seed: no.
- `scale`: 16 families, 3 generations, 3200 variants. Family-selective queries touch 18.3% of samples on seed 7. Gap on every seed: no.

A gap means PedPack and the no-pedigree ranker chose different layouts and
PedPack's measured fit objective is lower. The held-out byte ratio is
no-pedigree bytes divided by PedPack bytes on families that were not used
to choose the layout.

## Where the gap stops

On the three-generation pedigrees (`overlap`, cross rate 0.6) both rankers
select an all-samples layout. That is the exhaustive winner. Nuclear-family
duplication is 1.42, 1.33, and 1.43, and the family oracle is more expensive
than the variant oracle (6.344 vs 5.012, 4.732 vs 4.043, 6.612 vs 5.294).
Dropping the scan share from 0.35 to 0.20 (`results/stress/overlap_scan20`)
moves only seed 11 onto families (duplication 1.33, fit 3.316 vs 3.502).
Seeds 7 and 19 stay on the variant layout because the family layout is still
the more expensive one.

The larger cohort (`scale`, 102–147 samples, 3,200 variants) ties on seed 7,
where duplication of the family oracle is 1.45. Seeds 11 and 19 still pick
families, but the fit-objective ratios are only 1.018 and 1.001. Held-out
family bytes on those two seeds are 7.5× and 6.8× lower for the family layout.

The search is following the cost. The gap disappears when storing the
overlapping generation costs more, on the stated mix, than the selective
reads save.
