"""Fixed-length descriptors for the cost network.

Features are layout geometry, query shape, and a few summaries of the kinship
encoder. They do not include simulator outputs. Pedigree-specific columns can
be zeroed for the ablation that keeps the same candidate pool but hides
relatedness from the ranker.
"""

from __future__ import annotations

import numpy as np

from pedpack.cost import cover_group_indices
from pedpack.types import Cluster, Cohort, PreparedLayout, Query

FEATURE_NAMES: tuple[str, ...] = (
    "log_width",
    "log_mean_group",
    "log_n_groups",
    "group_size_cv",
    "trio_colocation",
    "pc_cut_rate",
    "mean_within_kinship",
    "duplication",
    "log_n_files",
    "log_n_blocks",
    "mean_file_over_block",
    "pack_naive",
    "order_group_major",
    "log_replication",
    "log_n_nodes",
    "q_scan_like",
    "q_log_sample_frac",
    "q_log_variant_frac",
    "q_log_n_trios",
    "q_trio_colocation",
    "q_sample_amp",
    "enc_within_cos",
    "enc_cut_cos",
    "enc_kept_cos",
    "pedigree_group",
)

PEDIGREE_FEATURES: tuple[str, ...] = (
    "trio_colocation",
    "pc_cut_rate",
    "mean_within_kinship",
    "q_trio_colocation",
    "enc_within_cos",
    "enc_cut_cos",
    "enc_kept_cos",
    "pedigree_group",
)


def pedigree_feature_index() -> list[int]:
    return [FEATURE_NAMES.index(name) for name in PEDIGREE_FEATURES]


def _cosine_members(embeddings: np.ndarray, members: list[int]) -> float:
    if len(members) < 2:
        return 0.0
    vecs = embeddings[members]
    norms = np.linalg.norm(vecs, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-8)
    unit = vecs / norms
    gram = unit @ unit.T
    mask = ~np.eye(len(members), dtype=bool)
    return float(gram[mask].mean())


def _encoder_summaries(layout: PreparedLayout, cohort: Cohort, embeddings: np.ndarray | None) -> tuple[float, float, float]:
    if embeddings is None:
        return 0.0, 0.0, 0.0
    within_num = 0.0
    within_den = 0.0
    for group in layout.groups:
        members = group.tolist()
        if len(members) < 2:
            continue
        within_num += _cosine_members(embeddings, members) * len(members)
        within_den += len(members)
    within = within_num / within_den if within_den else 0.0
    cut_vals: list[float] = []
    kept_vals: list[float] = []
    sets = layout.group_sets
    for child in range(cohort.n_samples):
        for parent in (int(cohort.father[child]), int(cohort.mother[child])):
            if parent < 0:
                continue
            cos = _pair_cosine(embeddings[parent], embeddings[child])
            if any(parent in s and child in s for s in sets):
                kept_vals.append(cos)
            else:
                cut_vals.append(cos)
    cut = float(np.mean(cut_vals)) if cut_vals else 0.0
    kept = float(np.mean(kept_vals)) if kept_vals else 0.0
    return within, cut, kept


def _pair_cosine(left: np.ndarray, right: np.ndarray) -> float:
    denom = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denom < 1e-12:
        return 0.0
    return float(np.dot(left, right) / denom)


def query_trio_colocation(layout: PreparedLayout, query: Query) -> float:
    if not query.trios:
        return 1.0
    hit = 0
    for father, mother, child in query.trios:
        if any(father in s and mother in s and child in s for s in layout.group_sets):
            hit += 1
    return hit / len(query.trios)


def feature_vector(
    layout: PreparedLayout,
    query: Query,
    cohort: Cohort,
    cluster: Cluster,
    embeddings: np.ndarray | None = None,
    encoder_stats: tuple[float, float, float] | None = None,
) -> np.ndarray:
    stats = layout.stats
    n_samples = max(1, cohort.n_samples)
    n_variants = max(1, cohort.n_variants)
    mean_file = stats.primary_bytes / max(1, stats.n_files)
    chosen, _uncovered = cover_group_indices(layout, query.samples)
    chosen_samples = 0
    for index in chosen:
        chosen_samples += int(layout.groups[index].size)
    amp = chosen_samples / max(1, len(query.samples))
    scan_like = 1.0 if len(query.samples) >= n_samples else 0.0
    if encoder_stats is None:
        enc_within, enc_cut, enc_kept = _encoder_summaries(layout, cohort, embeddings)
    else:
        enc_within, enc_cut, enc_kept = encoder_stats
    pedigree_group = 0.0 if stats.family in {"variant", "generic"} else 1.0
    values = [
        np.log1p(stats.width),
        np.log1p(stats.mean_group_size),
        np.log1p(stats.n_groups),
        stats.group_size_cv,
        stats.trio_colocation,
        stats.pc_cut_rate,
        stats.mean_within_kinship,
        stats.duplication,
        np.log1p(stats.n_files),
        np.log1p(stats.n_blocks),
        mean_file / cluster.block_size,
        1.0 if stats.pack == "naive" else 0.0,
        1.0 if stats.order == "group_major" else 0.0,
        np.log1p(cluster.replication),
        np.log1p(cluster.n_nodes),
        scan_like,
        np.log1p(len(query.samples) / n_samples),
        np.log1p(query.n_variants_span / n_variants),
        np.log1p(len(query.trios)),
        query_trio_colocation(layout, query),
        amp,
        enc_within,
        enc_cut,
        enc_kept,
        pedigree_group,
    ]
    vector = np.asarray(values, dtype=np.float64)
    if vector.shape != (len(FEATURE_NAMES),):
        raise RuntimeError("feature vector length drifted from FEATURE_NAMES")
    return vector


def feature_matrix(
    layouts: list[PreparedLayout],
    queries: list[Query],
    cohort: Cohort,
    cluster: Cluster,
    embeddings: np.ndarray | None = None,
    hide_pedigree: bool = False,
) -> np.ndarray:
    """Stack one row per (layout, query), layout-major.

    Encoder summaries do not depend on the query, so they are computed once
    per layout.
    """

    if not layouts or not queries:
        return np.zeros((0, len(FEATURE_NAMES)))
    rows: list[np.ndarray] = []
    for layout in layouts:
        enc_within, enc_cut, enc_kept = _encoder_summaries(layout, cohort, embeddings)
        for query in queries:
            rows.append(
                feature_vector(
                    layout,
                    query,
                    cohort,
                    cluster,
                    embeddings,
                    encoder_stats=(enc_within, enc_cut, enc_kept),
                )
            )
    matrix = np.vstack(rows)
    if hide_pedigree:
        matrix = matrix.copy()
        matrix[:, pedigree_feature_index()] = 0.0
    return matrix


def stack_targets(measurements: list) -> np.ndarray:
    """Four columns: log bytes, remote fraction, log shuffle, log milliseconds."""

    bytes_read = np.asarray([m.bytes_read for m in measurements], dtype=np.float64)
    remote = np.asarray([m.remote_fraction for m in measurements], dtype=np.float64)
    shuffle = np.asarray([m.shuffle_bytes for m in measurements], dtype=np.float64)
    millis = np.asarray([m.seconds * 1000.0 for m in measurements], dtype=np.float64)
    return np.column_stack(
        [
            np.log1p(np.maximum(bytes_read, 0.0)),
            np.clip(remote, 0.0, 1.0),
            np.log1p(np.maximum(shuffle, 0.0)),
            np.log1p(np.maximum(millis, 0.0)),
        ]
    )


def decode_targets(raw: np.ndarray) -> np.ndarray:
    """Invert ``stack_targets`` back to bytes, remote fraction, shuffle, seconds."""

    raw = np.asarray(raw, dtype=np.float64)
    decoded = np.zeros_like(raw)
    decoded[:, 0] = np.expm1(np.clip(raw[:, 0], 0.0, 40.0))
    decoded[:, 1] = np.clip(raw[:, 1], 0.0, 1.0)
    decoded[:, 2] = np.expm1(np.clip(raw[:, 2], 0.0, 40.0))
    decoded[:, 3] = np.expm1(np.clip(raw[:, 3], 0.0, 40.0)) / 1000.0
    return decoded
