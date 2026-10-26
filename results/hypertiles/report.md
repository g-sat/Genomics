# Hyper-tiles

`src/pedpack/hyper_tile.py`. Nine synthetic coordinates. Three pedigree groups are `(0,1,2)`, `(3,4,5)`, `(6,7,8)`. Cosine similarity at 0.5 links the same axis across families, so the phenotype components are `(0,3,6)`, `(1,4,7)`, `(2,5,8)`. Genotype payload is 100 bytes and phenotype payload is 20 bytes. A hyper-tile read pays 120 bytes per person in every opened component.

| Query | Pedigree tile | Phenotype tile | Split files |
| --- | ---: | ---: | ---: |
| Joint, one family | 360 | 1080 | 480 |
| Joint, one symptom component | 1080 | 360 | 960 |
| Genotype only, one family | 360 | 1080 | 300 |
| Phenotype only, one component | 1080 | 360 | 60 |

The Mendelian fraction is the share of queries that name a family. The rest name a symptom component.

| Workload | Fraction 0.8 | Fraction 0.5 | Fraction 0.2 |
| --- | --- | --- | --- |
| Both payloads | pedigree (504) | tie (720) | phenotype (504) |
| One payload | split (252) | split (180) | split (108) |

A joint family check is cheaper in a pedigree tile because the split layout opens three phenotype components to collect the second payload. A query that needs only one payload is cheaper in separate files, because a hyper-tile also reads the other payload. The coordinates are the test input. No clinical note, EHR row, or diagnostic platform was read.
