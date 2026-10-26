# Real chromosome 22 slice

2504 Phase 3 samples, 3000 biallelic sites, 605 pedigree families with two or more members.
The Phase 3 callset has no parent and child together. Family groups are the pedigree family ids.
Opaque reads are raw uint8 byte ranges rounded to 64 KB. Parquet is uncompressed column chunks.
The rule keeps a partition when stored bytes + 32 mean reads is below the all-samples tile.

Stored bytes: all-samples 7512000, nuclear 7512000 (1.000x), superpopulation 7512000 (1.000x).

| Workload | Mean opaque read (B) | Rule keeps the partition |
| --- | ---: | --- |
| pedigree family id, checks stay inside | 65536 | True |
| pedigree family id, check adds 12 other people | 851968 | False |
| pedigree family id, scans | 122683392 | False |
| superpopulation, checks stay inside | 188416 | True |
| superpopulation, checks name two populations | 365909 | True |
| superpopulation, scans | 7733248 | False |

All-samples window read: 524288 bytes.
One pedigree family, Parquet column chunks: 325 bytes.
Opaque all-samples / nuclear inside-check: 8.00x.
