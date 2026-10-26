# Chromosome 22 pointer stitch

`experiments/real_stitch.py`. High-coverage panel, first 3,000 biallelic SNVs, 3,202 samples, 602 trios. Dosages are counted GT alleles. Sample columns are packed in VCF order into 64 KB blocks. Reads use `FILE_FLAG_NO_BUFFERING`. Times are medians of 7 trials after one warmup read. The first nuclear family has three people, and those three columns fall in one block.

## 3,000 sites

| Read | Bytes moved | Median time |
| --- | ---: | ---: |
| Pointer stitch | 9000 logical, 12288 aligned | 0.108 ms |
| The one block that holds the trio | 65536 | 0.168 ms |
| Packed cohort file | 10027008 | 7.45 ms |

## 200 sites

| Read | Bytes moved | Median time |
| --- | ---: | ---: |
| Pointer stitch | 600 logical, 4096 aligned | 0.160 ms |
| The one block that holds the trio | 65536 | 0.331 ms |
| Packed cohort file | 655360 | 0.454 ms |

## Deletion alias

The withdrawn parent of that trio occupies one block. Replacing that block writes 65,536 bytes. Replacing the packed 3,000-site file writes 10,027,008 bytes. Median flush on this disk: 37.7 ms for the block, 22.2 ms for the file. Three people are each in two nuclear groups, of sizes 4 and 3. A family-tile deletion of one of them rewrites both groups.

## Parquet row groups

Uncompressed. One binary column, one row per sample.

| Layout | Bytes | Median `read` |
| --- | ---: | ---: |
| One row group, all samples | 9644279 column, file 9650656 | 10.5 ms |
| The trio's row group | 15069 | 0.755 ms |
| 2006 family and singleton row groups | file 34011188 |  |

The family file is larger because each row group has its own pages.
