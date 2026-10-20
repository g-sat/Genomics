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

## Same boundary at 600–1,200 samples

`configs/scale.json` keeps the locked cost model and the 5% family queries.
It draws 180 founder couples, 2 generations, and 8,000 variants. Seeds 7, 11,
and 19 produce 1,207, 1,175, and 1,176 people. Nuclear duplication is
1.28–1.30. PedPack selects nuclear families and matches the family oracle.
The no-pedigree ranker selects all-samples. Fit objectives without pedigree
features are 1.20×, 1.22×, and 1.28× PedPack. Held-out family bytes are
16.3×, 16.7×, and 17.6×. Scan bytes go the other way, 1.15–1.17×, because
the family layout stores the overlapping parents twice.

The fit objective already contains build, storage, and query time. On these
three seeds the family layout costs more to build and more to store, and
less to query:

| Seed | People | PedPack build / storage / query (s) | No-pedigree build / storage / query (s) |
| ---: | ---: | --- | --- |
| 7 | 1207 | 79.25 / 4.90 / 5.72 | 71.41 / 4.26 / 7.81 |
| 11 | 1175 | 76.77 / 4.74 / 5.89 | 69.51 / 4.15 / 7.81 |
| 19 | 1176 | 76.32 / 4.70 / 4.35 | 69.57 / 4.15 / 6.53 |

`configs/scale_overlap.json` is the other side of the boundary at the same
variant count: 80 founder couples, 3 generations, cross rate 0.6. The draws
have 667, 675, and 685 people. Family-oracle duplication is 1.42–1.44, and
that oracle is more expensive than all-samples (seed 7: 218.7 vs 173.1).
Both rankers select all-samples. The gap closes for the same reason it
closed on the 70–100 person pedigrees.

## Same open side at about 3,300 samples

`configs/scale_3k.json` uses 500 founder couples, 2 generations, cross rate
0.25, and 8,000 variants. Seeds 7, 11, and 19 have 3,370, 3,361, and 3,328
people. Nuclear duplication is 1.293–1.298. PedPack selects nuclear families
and matches the family oracle. The no-pedigree ranker selects all-samples.
Fit objectives without pedigree features are 1.244×, 1.211×, and 1.208×
PedPack. Held-out family bytes are 16.8×, 17.3×, and 16.7×. Scan bytes are
1.16–1.17× higher for the family layout. Build and storage are higher for
families; query time is lower. The split did not move between 1,200 and
3,300 people.

## 1000 Genomes pedigree, simulated dosages

`data/pedigrees/1kg_complete_parents.ped` is the phase 3 pedigree after
dropping people with only one recorded parent. It has 3,375 people, 636
children, and 2,739 founders. `configs/public_1kg.json` runs the same
selector and the same 8,000 simulated variants. Nuclear duplication is
1.001. PedPack and the blind ranker both reject all-samples. Their fit
objectives agree to within 0.2% on every seed. The all-samples fit
objective is 1.42–1.49× PedPack, and held-out family bytes are 18.6–20.1×
higher. On a pedigree that does not overlap, the group-size term is enough.

## The two-layout surface

`family_versus_variant` is that trade-off with the search removed. At 1,200
samples and 8,000 variants the all-samples objective divided by the family
objective is 1.06 at duplication 1.30 and a 5% family query, 0.97 at
duplication 1.43 and the same query, and 0.83 at duplication 1.30 and a 60%
family query. Family tiles win on 32 of the 80 cells. The picture is
`results/boundary/heatmap.png`.

Simulator checks of only those two layouts, at 1,600 variants, agree on the
sign in 8 of 12 draws. The four misses are all the same direction: the
formula has already switched to all-samples, and the simulator still picks
nuclear families, by 12% at one 15% query and by 0.3–4.7% at 40%. The
formula is early on that edge. It does not invent a family win the
simulator rejects.

## Hail, one virtual machine

`results/hail/report.md` is the first local-mode run. The family object
there is an exploded Table, and it is about 20× the MatrixTable even though
duplication is 1.30. That comparison prices the encoding.

`results/hail_columnar/report.md` is the fairer one, seed 7 only. Both sides
are integer MatrixTables. The all-samples table has 16 locus partitions,
builds in 18.7 s, and occupies 3.0 MB. The family side is 361 MatrixTables,
builds in 937 s, and occupies 33.9 MB. A held-out query that returns 250,518
entries takes 0.35–1.7 s on the single table and 7.1–8.1 s across the 18
covering family tables. The scan takes 0.55 s against 169 s. The simulator
still prefers the family layout on bytes, because a file open is 1 ms there.
Hail's cost at this size is the Spark job, and the many small tables lose.

The first run's numbers, kept so the encoding effect stays visible:

| Seed | MatrixTable build s | Family build s | MatrixTable bytes | Family bytes | Held-out MatrixTable s | Held-out family s |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 12.8 | 131.5 | 2.96e6 | 5.80e7 | 0.718 | 0.359 |
| 11 | 4.6 | 118.8 | 2.88e6 | 5.60e7 | 0.429 | 0.348 |
| 19 | 4.4 | 123.9 | 2.89e6 | 5.57e7 | 0.355 | 0.361 |

The family table is about 20× larger than the MatrixTable. Nuclear
duplication is only 1.28–1.30, so the extra bytes are the exploded row
encoding, not the extra parent copies. Held-out query means favor the
family key on seeds 7 and 11 and tie on seed 19. Individual queries cross,
and every query finishes in under a second. The scan is not consistently
slower for the family table (seed 11: 0.66 s against 0.89 s). At this size
Hail's MatrixTable is the better full-cost object, and the simulator's
16× selective-byte gap is not visible in the Spark timings.

## Packed Parquet, the layout the search chose

`results/parquet/report.md` is seed 7 again, Zstandard, 32 groups per file,
12 files. That is `count:32`, not one object per family. Build is 1.69 s /
4.54 MB for all-samples and 1.28 s / 5.99 MB for the packed files. Held-out
compressed bytes are 84.3 KB against 85.2 KB. Packed queries are slower
because they open 8–10 files. The scan is 1.30× the bytes. Column projection
on the single file already skips the people the query did not name, which is
the saving the opaque-tile simulator was charging the all-samples layout for.

## HDFS, one host

`results/hdfs/report.md` is a Hadoop 3.2.1 cluster of one NameNode and three
DataNodes on this machine, block size 64 KiB, replication 3. On select_05 the
client received 14–20× fewer held-out bytes from the family layout than from
the all-samples layout, and the NameNode mapped those reads onto about 6–8
blocks instead of 75–93. On the three-generation overlap regime the two
layouts are the same all-samples files, and the bytes match. The client is
not a DataNode, so this run has no remote-read fraction. It is not three
separate machines. The compose project was stopped after that measurement.
Hail is the local-mode run above, not a job on this Hadoop cluster.
