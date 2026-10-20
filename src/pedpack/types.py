"""Shared records for cohorts, queries, layouts, and measurements."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class Cluster:
    """Block-storage parameters used by the cost model and the simulator.

    Bandwidths are bytes per second. ``storage_scan_equivalents`` converts
    replicated bytes into the same objective as query time: storing the
    reference layout is priced as this many reference full scans.
    """

    n_nodes: int = 8
    block_size: int = 262_144
    replication: int = 3
    bytes_per_genotype: float = 48.0
    compress_base: float = 0.45
    kinship_alpha: float = 0.55
    bw_local: float = 2.0e8
    bw_remote: float = 8.0e7
    bw_shuffle: float = 1.0e8
    t_open: float = 0.004
    genotypes_per_second: float = 2.0e7
    failed_nodes: tuple[int, ...] = ()
    storage_scan_equivalents: float = 0.25

    def with_updates(self, **kwargs: object) -> Cluster:
        data = {
            "n_nodes": self.n_nodes,
            "block_size": self.block_size,
            "replication": self.replication,
            "bytes_per_genotype": self.bytes_per_genotype,
            "compress_base": self.compress_base,
            "kinship_alpha": self.kinship_alpha,
            "bw_local": self.bw_local,
            "bw_remote": self.bw_remote,
            "bw_shuffle": self.bw_shuffle,
            "t_open": self.t_open,
            "genotypes_per_second": self.genotypes_per_second,
            "failed_nodes": self.failed_nodes,
            "storage_scan_equivalents": self.storage_scan_equivalents,
        }
        data.update(kwargs)
        return Cluster(**data)  # type: ignore[arg-type]


@dataclass
class Cohort:
    """Genotype matrix aligned to a pedigree.

    Column ``j`` is one individual. ``father[j]`` and ``mother[j]`` are column
    indices of the parents, or -1 for a founder. ``kinship`` is the n x n
    matrix of kinship coefficients. Columns are stored in a permuted order so
    that consecutive columns are not nuclear families.
    """

    ids: np.ndarray
    positions: np.ndarray
    genotypes: np.ndarray
    father: np.ndarray
    mother: np.ndarray
    sex: np.ndarray
    family_id: np.ndarray
    kinship: np.ndarray

    @property
    def n_samples(self) -> int:
        return int(self.genotypes.shape[1])

    @property
    def n_variants(self) -> int:
        return int(self.genotypes.shape[0])


@dataclass(frozen=True)
class Query:
    """One hereditary-analysis request against the matrix.

    ``samples`` and ``trios`` are column indices. A trio is
    ``(father, mother, child)``. ``workload`` is the protocol bucket
    (``fit``, ``scan``, ``family_selective``, ``extended``, ``heldout_family``,
    ``heldout_extended``).
    """

    name: str
    workload: str
    samples: tuple[int, ...]
    v0: int
    v1: int
    trios: tuple[tuple[int, int, int], ...]

    @property
    def n_variants_span(self) -> int:
        return max(0, self.v1 - self.v0)


@dataclass(frozen=True)
class Candidate:
    """One point in the layout space, before bytes are assigned."""

    label: str
    family: str
    groups: tuple[tuple[int, ...], ...]
    width: int
    pack: str
    order: str


@dataclass
class LayoutStats:
    """Geometry the analytical model is allowed to see."""

    n_groups: int
    n_tiles: int
    n_files: int
    n_blocks: int
    storage_bytes: float
    primary_bytes: float
    mean_group_size: float
    group_size_cv: float
    trio_colocation: float
    pc_cut_rate: float
    mean_within_kinship: float
    duplication: float
    width: int
    order: str
    pack: str
    family: str
    build_seconds: float
    label: str


@dataclass
class PreparedLayout:
    """Tile byte ranges inside packed files. Independent of node placement."""

    candidate: Candidate
    groups: list[np.ndarray]
    group_sets: list[set[int]]
    tile_group: np.ndarray
    tile_v0: np.ndarray
    tile_v1: np.ndarray
    tile_file: np.ndarray
    tile_offset: np.ndarray
    tile_nbytes: np.ndarray
    file_nbytes: np.ndarray
    stats: LayoutStats
    group_kinship: np.ndarray = field(default_factory=lambda: np.zeros(0))


@dataclass(frozen=True)
class Measurement:
    bytes_read: float
    remote_bytes: float
    shuffle_bytes: float
    seconds: float
    files_touched: int
    blocks_touched: int
    uncovered: int
    unavailable_blocks: int

    @property
    def remote_fraction(self) -> float:
        if self.bytes_read <= 0:
            return 0.0
        return self.remote_bytes / self.bytes_read
