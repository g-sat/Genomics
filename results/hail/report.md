# Hail baseline, local mode

Hail 0.2.139 in Spark local mode on one WSL virtual machine. This is
not a multi-node Hail or HDFS cluster. The default object is a
variant-keyed MatrixTable. The family object is a Hail Table keyed by
nuclear family, so a person who sits in two families is stored twice.
Both reads use interval filters. Build time is the Hail write. Storage
is the directory size. Query time is the Spark action for that read.
On these matrices the MatrixTable is one partition, so a variant window
still opens the whole file. The family table is 64 partitions.

| Seed | Samples | Variants | Duplication | MT build s | Family build s | MT bytes | Family bytes | Held-out MT s | Held-out family s |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 1207 | 8000 | 1.300 | 12.8 | 131.5 | 2962103 | 57985834 | 0.718 | 0.359 |
| 11 | 1175 | 8000 | 1.291 | 4.6 | 118.8 | 2882511 | 55997500 | 0.429 | 0.348 |
| 19 | 1176 | 8000 | 1.281 | 4.4 | 123.9 | 2887379 | 55702898 | 0.355 | 0.361 |

## Queries

### Seed 7

| Query | Samples | Variants | Groups | MT s | MT entries | Family s | Family entries |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| hold-fam-0 | 86 | 2913 | 18 | 0.939 | 250518 | 0.324 | 250518 |
| hold-fam-1 | 71 | 2907 | 18 | 0.669 | 206397 | 0.389 | 209304 |
| hold-fam-2 | 74 | 2103 | 18 | 0.854 | 155622 | 0.392 | 157725 |
| hold-fam-3 | 78 | 1677 | 18 | 0.410 | 130806 | 0.331 | 132483 |
| scan-0 | 1207 | 4425 | 361 | 0.546 | 5340975 | 0.794 | 6942825 |

### Seed 11

| Query | Samples | Variants | Groups | MT s | MT entries | Family s | Family entries |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| hold-fam-0 | 76 | 1952 | 18 | 0.408 | 148352 | 0.313 | 150304 |
| hold-fam-1 | 80 | 2134 | 18 | 0.388 | 170720 | 0.360 | 172854 |
| hold-fam-2 | 76 | 1710 | 18 | 0.497 | 129960 | 0.343 | 131670 |
| hold-fam-3 | 76 | 2253 | 18 | 0.423 | 171228 | 0.376 | 171228 |
| scan-0 | 1175 | 6614 | 351 | 0.894 | 7771450 | 0.656 | 10033438 |

### Seed 19

| Query | Samples | Variants | Groups | MT s | MT entries | Family s | Family entries |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| hold-fam-0 | 71 | 2427 | 17 | 0.327 | 172317 | 0.645 | 177171 |
| hold-fam-1 | 75 | 2226 | 17 | 0.359 | 166950 | 0.229 | 166950 |
| hold-fam-2 | 69 | 2965 | 17 | 0.285 | 204585 | 0.384 | 204585 |
| hold-fam-3 | 74 | 1727 | 17 | 0.446 | 127798 | 0.185 | 127798 |
| scan-0 | 1176 | 6452 | 345 | 0.420 | 7587552 | 0.888 | 9716712 |
