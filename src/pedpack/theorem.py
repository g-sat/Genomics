"""When a family layout has the lower objective, in the large-payload limit.

The simulator and ``family_versus_variant`` keep block rounding and file
opens. This module drops both, which is valid once every read spans many
blocks. The resulting inequality is closed form. On the 1,200 by 8,000
grid used in the paper it has the same sign as ``family_versus_variant``
on 190 of 192 cells. The two misses sit within half a percent of a tie.

Assumptions, all required for the inequality:

- Both layouts are charged the same remote fraction.
- Kinship does not change the stored size.
- Queried trios sit inside the groups that are opened, so shuffle is zero.
- A family query reads a fraction ``f`` of the stored family bytes, and an
  extended query reads ``min(1, 2f)``.
- The family layout stores ``D`` bytes per all-samples byte.
- Build reads the source matrix once either way. Only the write of the
  primary copy and its replicas differs.
- Storage of the all-samples layout is priced at ``storage_scan_equivalents``
  full scans of that layout.

Let ``S`` be the all-samples payload and ``c_byte`` the seconds per payload
byte at the shared remote fraction. Scan time tends to ``c_byte * S`` and
``c_byte * D * S``. A window that covers a fraction ``w`` of the variants
scales the selective reads by ``w``. Write cost per extra payload byte is
``c_write``. Dividing the objective gap by ``S > 0`` leaves

    (D - 1) * (c_write + equivalents * c_byte) + horizon * c_byte * Q(D, f)

with

    Q = scan*(D - 1) + w*family*(f*D - 1) + w*extended*(min(1, 2f)*D - 1).

The family layout is cheaper when that quantity is negative.
"""

from __future__ import annotations

import math

from pedpack.types import Cluster


def byte_rate(cluster: Cluster, remote_fraction: float) -> float:
    """Seconds per payload byte at a fixed remote fraction, excluding opens."""

    remote = min(1.0, max(0.0, float(remote_fraction)))
    return (
        (1.0 - remote) / cluster.bw_local
        + remote / cluster.bw_remote
        + 1.0 / (cluster.bytes_per_genotype * cluster.genotypes_per_second)
    )


def write_rate(cluster: Cluster) -> float:
    """Seconds to write one extra payload byte, primary local and replicas remote."""

    replicas = max(0, cluster.replication - 1)
    return 1.0 / cluster.bw_local + replicas / cluster.bw_remote


def query_mix(
    duplication: float,
    family_fraction: float,
    window_fraction: float = 0.30,
    scan_share: float = 0.35,
    family_share: float = 0.40,
    extended_share: float = 0.25,
) -> float:
    """Query-time gap per all-samples scan byte. Negative means family is cheaper to read."""

    extended = min(1.0, 2.0 * family_fraction)
    return (
        scan_share * (duplication - 1.0)
        + window_fraction * family_share * (family_fraction * duplication - 1.0)
        + window_fraction * extended_share * (extended * duplication - 1.0)
    )


def objective_gap(
    duplication: float,
    family_fraction: float,
    cluster: Cluster,
    remote_fraction: float = 0.60,
    window_fraction: float = 0.30,
    horizon: float = 32.0,
    scan_share: float = 0.35,
    family_share: float = 0.40,
    extended_share: float = 0.25,
) -> float:
    """Family objective minus all-samples objective, per all-samples payload byte.

    Negative means the family layout wins in the large-payload limit.
    """

    if duplication < 1:
        raise ValueError("duplication is at least 1")
    per_byte = byte_rate(cluster, remote_fraction)
    mix = query_mix(
        duplication,
        family_fraction,
        window_fraction=window_fraction,
        scan_share=scan_share,
        family_share=family_share,
        extended_share=extended_share,
    )
    storage = cluster.storage_scan_equivalents * per_byte
    return (duplication - 1.0) * (write_rate(cluster) + storage) + horizon * per_byte * mix


def block_ratio(numerator: int, denominator: int, block: int) -> tuple[int, int, float]:
    """Blocks read for two payloads, and the ratio of those block counts.

    A payload of zero blocks is reported as zero. The ratio is undefined then
    and returned as 0.
    """

    if block < 1:
        raise ValueError("block size must be positive")
    if numerator < 0 or denominator < 0:
        raise ValueError("payloads must be non-negative")
    num_blocks = math.ceil(numerator / block) if numerator else 0
    den_blocks = math.ceil(denominator / block) if denominator else 0
    if den_blocks == 0:
        return num_blocks, den_blocks, 0.0
    return num_blocks, den_blocks, num_blocks / den_blocks
