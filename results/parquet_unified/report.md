# One Parquet file for the packed families

Seed 7 of the scaled cohort, 1,207 people, 361 nuclear groups, duplication
1.300. The all-samples file and the 12-file pack are the same sizes as
`results/parquet/` (4,538,050 and 5,993,701 bytes). The queries in this
file are this process's held-out draw. They are not the four windows in
`results/parquet/report.md`, so the two reports are not row-aligned.
Compression is Zstandard.
All-samples is one file, one column per person, row groups of
800 variants. Split is the search pack:
32 nuclear groups per file,
12 files, duplication 1.300.
Unified is those same duplicated columns in one file, same row-group
width, so a query opens one footer and still skips row groups outside
its variant window. Wide is that one file with a single row group
covering every variant, so a window cannot skip.

| Layout | Files | Row group | Build s | Bytes |
| --- | ---: | ---: | ---: | ---: |
| all-samples | 1 | 800 | 9.61 | 4538050 |
| split count:32 | 12 | 800 | 2.88 | 5993701 |
| unified | 1 | 800 | 2.67 | 6013677 |
| unified wide | 1 | 8000 | 1.00 | 2670838 |

| Query | Groups | Split files | All-samples s | All-samples B | Split s | Split B | Unified s | Unified B | Wide s | Wide B |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| hold-fam-0 | 23 | 11 | 0.096 | 109545 | 0.217 | 110630 | 0.059 | 110630 | 0.036 | 151336 |
| hold-fam-1 | 33 | 11 | 0.043 | 151753 | 0.066 | 153941 | 0.063 | 153941 | 0.021 | 211781 |
| hold-fam-2 | 46 | 11 | 0.042 | 204966 | 0.071 | 209393 | 0.061 | 209393 | 0.029 | 286613 |
| hold-fam-3 | 45 | 12 | 0.050 | 258337 | 0.069 | 269270 | 0.068 | 269270 | 0.040 | 295818 |
| scan-0 | 361 | 12 | 0.149 | 2310991 | 0.194 | 3004157 | 0.170 | 3004157 | 0.069 | 2355559 |

Mean held-out: all-samples 0.058 s / 181150 B; split 0.106 s / 185808 B; unified 0.063 s / 185808 B; wide 0.031 s / 236387 B.
