# Hail baseline, local mode

Hail 0.2.139 in Spark local mode on one WSL virtual machine. This is
not a multi-node Hail or HDFS cluster. Both objects are MatrixTables
with the same integer entry type. The default table is variant-keyed
and repartitioned along the locus axis. The family side is one
MatrixTable per nuclear family, so a person in two families is stored
twice and a family query opens only the covering tables.

| Seed | Samples | Variants | Duplication | MT build s | Family build s | MT bytes | Family bytes | Held-out MT s | Held-out family s |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 1207 | 8000 | 1.300 | 18.7 | 936.6 | 2954268 | 33884477 | 0.763 | 7.596 |

## Queries

### Seed 7

| Query | Samples | Variants | Groups | MT s | MT entries | Family s | Family entries |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| hold-fam-0 | 86 | 2913 | 18 | 1.714 | 250518 | 8.053 | 250518 |
| hold-fam-1 | 71 | 2907 | 18 | 0.584 | 206397 | 7.434 | 209304 |
| hold-fam-2 | 74 | 2103 | 18 | 0.348 | 155622 | 7.842 | 157725 |
| hold-fam-3 | 78 | 1677 | 18 | 0.408 | 130806 | 7.055 | 132483 |
| scan-0 | 1207 | 4425 | 361 | 0.552 | 5340975 | 169.162 | 6942825 |
