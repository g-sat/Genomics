"""Encoder and cost-network training loops."""

from __future__ import annotations

import numpy as np

from pedpack.genotype import build_cohort
from pedpack.models import train_embeddings, train_surrogate
from pedpack.pedigree import Pedigree


def _two_families() -> "object":
    pedigree = Pedigree(
        ids=np.array(["a", "b", "c", "d", "e", "f"], dtype=object),
        father=np.array([-1, -1, 0, -1, -1, 3], dtype=np.int32),
        mother=np.array([-1, -1, 1, -1, -1, 4], dtype=np.int32),
        sex=np.array([1, 2, 1, 1, 2, 2], dtype=np.int32),
        family_id=np.array([0, 0, 0, 1, 1, 1], dtype=np.int32),
    )
    return build_cohort(pedigree, n_variants=30, rng=np.random.default_rng(0), permute_columns=False)


def test_encoder_puts_parent_and_child_closer_than_unrelated_founders():
    cohort = _two_families()
    embedded = train_embeddings(cohort, steps=200, hidden=16, seed=0)

    def cosine(i: int, j: int) -> float:
        denom = float(np.linalg.norm(embedded[i]) * np.linalg.norm(embedded[j]))
        return float(np.dot(embedded[i], embedded[j]) / denom)

    parent_child = cosine(0, 2)
    unrelated = cosine(0, 3)
    assert parent_child > unrelated + 0.05


def test_cost_net_fits_a_smooth_target():
    rng = np.random.default_rng(0)
    features = rng.normal(size=(100, 5))
    direction = np.array([0.4, -0.2, 0.15, 0.0, 0.25])
    signal = features @ direction
    targets = np.column_stack(
        [
            signal,
            np.clip(0.5 + 0.05 * features[:, 0], 0.0, 1.0),
            signal + 0.3,
            signal - 0.2,
        ]
    )
    model = train_surrogate(features, targets, hidden=16, epochs=80, seed=0, patience=20)
    predicted = model.predict_raw(features)
    correlation = float(np.corrcoef(predicted[:, 0], targets[:, 0])[0, 1])
    assert correlation > 0.85
