# Where the pedigree gap holds

The core positive results are three re-fits of the locked generator. Each one
is a complete config, not an overlay someone has to reconstruct.

| Config | What changes from `configs/pedigree_gap.json` | Result on seeds 7, 11, 19 |
| --- | --- | --- |
| `configs/positive/select_05.json` | Family queries touch 5% of nuclear families | PedPack picks nuclear families. The no-pedigree ranker picks all-samples. Fit objective 1.12–1.19× higher without pedigree features. Held-out family bytes 12–16× higher. |
| `configs/positive/select_15.json` | Family queries touch 15% of nuclear families | Same split of layouts. Fit objective 1.11–1.21×. Held-out bytes 5.4–6.0×. |
| `configs/positive/mix_50.json` | Fit mix is 50% scans and 50% family queries. No extended queries in the fit. Family fraction stays the locked 5–15% draw. | Same split of layouts. Fit objective 1.04–1.10×. Held-out bytes 5.0–8.8×. |

Shared generator, unless a cell above says otherwise: 10 founder couples, 2
generations, 2.2 children per couple, cross rate 0.25, at least 48 people,
1,600 variants. Cluster: 8 nodes in the simulator, 64 KB blocks, replication
3, 400 bytes per genotype. Horizon 32. Storage priced at 0.25 reference
scans. The numbers are in `results/stress/select_05`, `select_15`, and
`mix_50`. The figure is `results/stress/boundary.png`.

## Where it closes

| Config | Generator | What the search does |
| --- | --- | --- |
| `results/stress/select_40` | Same pedigrees. Queries touch 40% of families (about half the samples, because families overlap). | Layouts still differ. Fit ratio falls to 1.00–1.13. Seed 19 is 3.686 vs 3.701. Held-out bytes are still 2.3–2.8×. |
| `results/stress/overlap` | 3 generations, cross rate 0.6. Nuclear duplication 1.33–1.43. | Both rankers select all-samples. The family oracle is more expensive (for example 6.344 vs 5.012 on seed 7). |
| `results/stress/overlap_scan20` | Same 3-generation pedigrees, scan share lowered to 0.20. | Seed 11, duplication 1.33, moves to families (3.316 vs 3.502). Seeds 7 and 19 stay on all-samples because the family layout still costs more. |
| `results/stress/scale` | 16 founder couples, 3 generations, cross rate 0.45, 3,200 variants, 102–147 people. | Seed 7, family-oracle duplication 1.45, selects all-samples. Seeds 11 and 19 select families with fit ratios 1.018 and 1.001. |

The boundary is duplication against the scan share, not a failure to see the
pedigree. Around duplication 1.21–1.27 and a scan share of 0.35–0.50, family
tiles win the selective queries by enough to pay for the extra copies. Around
duplication 1.4, the extra copies cost more than those queries save, and both
rankers leave the family layout.

## HDFS, one host

`results/hdfs/report.md` is a Hadoop 3.2.1 cluster of one NameNode and three
DataNodes on this machine, block size 64 KiB, replication 3. On select_05 the
client received 14–20× fewer held-out bytes from the family layout than from
the all-samples layout, and the NameNode mapped those reads onto about 6–8
blocks instead of 75–93. On the three-generation overlap regime the two
layouts are the same all-samples files, and the bytes match. The client is
not a DataNode, so this run has no remote-read fraction. It is not three
separate machines, and it does not include Hail.
