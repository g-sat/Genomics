# High-coverage chromosome 22, real trios

3202 samples, 3000 biallelic SNVs from 3290 variant rows.
602 parent-child pairs have both parents in the VCF.
597 nuclear families have a sequenced child. Stored slots / samples = 1.001.
The inside check is the first such family (3 people) over 200 sites.

| Block | All-samples window | Trio window | Trio+12 | All-samples scan | Family-file scan | Window ratio |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 65536 | 655360 | 65536 | 851968 | 9633792 | 131465216 | 10.00 |
| 134217728 | 134217728 | 134217728 | 1744830464 | 134217728 | 269240762368 | 1.00 |

The 64 KB trio read is one block. The all-samples window is ten blocks. Adding twelve outsiders opens thirteen blocks, 851,968 bytes, which is more than the all-samples window, so that check does not keep the family files.

At 128 MB the trio and the all-samples window are each one block. The trio-plus-12 and family-file scan counts are one full block per opened file. They are the cost of leaving each family in its own file, not extra genotype bytes.

Duplication is 1.001. The 602 trios are almost disjoint, so this callset does not exercise the overlapping-parent penalty.
