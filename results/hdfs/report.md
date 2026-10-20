# HDFS measurements

One Windows host, Docker, Hadoop 3.2.1. One NameNode and three DataNodes
(`deploy/hdfs-compose.yml`). Replication 3. Block size 64 KiB, which required
lowering `dfs.namenode.fs-limits.min-block-size` from its 1 MiB default.
The client was another container on the same Docker network. It is not a
DataNode, so every block was fetched from a DataNode. That is not the
simulator's remote-read fraction, and this file does not report one.

Files are the layouts chosen on the simulator: nuclear-family layouts where
PedPack picked them, the all-samples layout from the no-pedigree ranker, and
the generic oracle. Each genotype is padded to the modeled tile size. There
is no codec. Queries are the same family-selective and held-out family lists
used in the simulator. Bytes below are what WebHDFS returned to the client.
Blocks are the NameNode's block locations for those ranges. Storage is
`spaceConsumed`, which includes the three replicas.

## Held-out family reads

| Scenario | Seed | PedPack bytes | No-pedigree bytes | Byte ratio | PedPack seconds | No-pedigree seconds | PedPack blocks | No-pedigree blocks |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| select_05 | 7 | 3.51e5 | 5.76e6 | 16.4 | 0.012 | 0.077 | 7.7 | 93.2 |
| select_05 | 11 | 3.19e5 | 4.53e6 | 14.2 | 0.022 | 0.103 | 5.8 | 74.7 |
| select_05 | 19 | 2.68e5 | 5.42e6 | 20.2 | 0.019 | 0.082 | 5.7 | 88.5 |
| select_15 | 7 | 9.8e5 | 6.81e6 | 6.9 | 0.038 | 0.083 | 23.2 | 110.2 |
| select_15 | 11 | 8.1e5 | 5.28e6 | 6.5 | 0.015 | 0.070 | 13.7 | 86.8 |
| select_15 | 19 | 8.5e5 | 5.73e6 | 6.7 | 0.021 | 0.071 | 16.5 | 93.7 |
| mix_50 | 7 | 8.6e5 | 6.77e6 | 7.9 | 0.027 | 0.080 | 17.8 | 109.3 |
| mix_50 | 11 | 9.1e5 | 5.29e6 | 5.8 | 0.018 | 0.071 | 16.3 | 86.7 |
| mix_50 | 19 | 5.8e5 | 6.10e6 | 10.5 | 0.017 | 0.076 | 12.3 | 99.3 |
| overlap | 7 | 7.66e6 | 7.66e6 | 1.0 | 0.062 | 0.073 | 118.5 | 121.7 |
| overlap | 11 | 6.77e6 | 6.77e6 | 1.0 | 0.079 | 0.072 | 109.8 | 109.8 |
| overlap | 19 | 8.35e6 | 8.35e6 | 1.0 | 0.103 | 0.119 | 133.0 | 133.0 |

On the two-generation selective regimes the family layout returns fewer bytes,
opens fewer blocks, and finishes the client read sooner. On the
three-generation overlap regime both rankers chose an all-samples layout, and
the NameNode serves the same bytes.

Replicated storage is slightly higher for the family layouts. On select_05 the
three seeds use 72.3, 53.2, and 64.9 MB with the family layout, against 64.6,
49.4, and 56.9 MB for the all-samples layout. File counts for the family
layouts are 5, 1, and 3. The all-samples layouts are one file each.

The generic layout also stores every sample once. Its held-out bytes are in
the same range as the all-samples layout. The simulator charges that layout
for shuffle when trios are split. This client does not run that shuffle, so
the generic penalty does not appear here.

The cohorts are the same 58–147 people as the simulator. This is not a
multi-machine cluster and it is not a Hail comparison.
