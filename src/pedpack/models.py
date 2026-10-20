"""Kinship encoder and the query-cost network.

The encoder is trained only on the pedigree: cosine similarity is regressed
onto scaled kinship. It never sees query times.

The cost network is a small multilayer perceptron. It maps layout and query
features to the simulator's four execution measurements. Storage and build
stay in the closed-form model, because those quantities are known once the
files have been described.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
from torch import nn

from pedpack.features import decode_targets
from pedpack.types import Cohort


def node_features(cohort: Cohort) -> np.ndarray:
    n = cohort.n_samples
    children = np.zeros(n, dtype=np.float64)
    for child in range(n):
        if int(cohort.father[child]) >= 0:
            children[int(cohort.father[child])] += 1.0
            children[int(cohort.mother[child])] += 1.0
    founders = ((cohort.father < 0) & (cohort.mother < 0)).astype(np.float64)
    n_parents = (cohort.father >= 0).astype(np.float64) + (cohort.mother >= 0).astype(np.float64)
    return np.column_stack(
        [
            founders,
            cohort.sex.astype(np.float64) / 2.0,
            np.log1p(children),
            n_parents,
            cohort.kinship.mean(axis=1),
        ]
    ).astype(np.float32)


def pedigree_edges(cohort: Cohort) -> torch.Tensor:
    """Undirected parent-child and mate edges, stored as a 2 x E index."""

    src: list[int] = []
    dst: list[int] = []
    mates: set[tuple[int, int]] = set()
    for child in range(cohort.n_samples):
        parents = [int(cohort.father[child]), int(cohort.mother[child])]
        live = [p for p in parents if p >= 0]
        for parent in live:
            src.extend((parent, child))
            dst.extend((child, parent))
        if len(live) == 2:
            a, b = (live[0], live[1]) if live[0] < live[1] else (live[1], live[0])
            mates.add((a, b))
    for a, b in mates:
        src.extend((a, b))
        dst.extend((b, a))
    if not src:
        return torch.zeros((2, 0), dtype=torch.long)
    return torch.tensor([src, dst], dtype=torch.long)


def _mean_aggregate(hidden: torch.Tensor, edge_index: torch.Tensor) -> torch.Tensor:
    n, _dim = hidden.shape
    if edge_index.numel() == 0:
        return torch.zeros_like(hidden)
    src, dst = edge_index
    out = torch.zeros_like(hidden)
    deg = torch.zeros(n, 1, dtype=hidden.dtype)
    out.index_add_(0, dst, hidden[src])
    ones = torch.ones(src.shape[0], 1, dtype=hidden.dtype)
    deg.index_add_(0, dst, ones)
    return out / deg.clamp(min=1.0)


class KinshipEncoder(nn.Module):
    """Two rounds of mean aggregation, then a linear map to the embedding."""

    def __init__(self, in_dim: int, hidden: int) -> None:
        super().__init__()
        self.lift = nn.Linear(in_dim, hidden)
        self.message = nn.Linear(hidden * 2, hidden)
        self.project = nn.Linear(hidden, hidden)

    def forward(self, features: torch.Tensor, edge_index: torch.Tensor) -> torch.Tensor:
        hidden = torch.relu(self.lift(features))
        for _ in range(2):
            pooled = _mean_aggregate(hidden, edge_index)
            hidden = torch.relu(self.message(torch.cat([hidden, pooled], dim=1)))
        return self.project(hidden)


def train_embeddings(
    cohort: Cohort,
    steps: int = 120,
    hidden: int = 32,
    seed: int = 0,
    lr: float = 1e-2,
) -> np.ndarray:
    """Return one L2-normalized embedding per sample."""

    torch.manual_seed(seed)
    features = torch.tensor(node_features(cohort))
    edges = pedigree_edges(cohort)
    target_np = np.clip(cohort.kinship / 0.5, 0.0, 1.0).astype(np.float32)
    target = torch.tensor(target_np)
    model = KinshipEncoder(features.shape[1], hidden)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    model.train()
    for _ in range(steps):
        optimizer.zero_grad()
        embedded = torch.nn.functional.normalize(model(features, edges), dim=1, eps=1e-8)
        cosine = embedded @ embedded.T
        loss = torch.mean((cosine - target) ** 2)
        loss.backward()
        optimizer.step()
    model.eval()
    with torch.no_grad():
        embedded = torch.nn.functional.normalize(model(features, edges), dim=1, eps=1e-8)
    return embedded.numpy().astype(np.float64)


class CostNet(nn.Module):
    def __init__(self, in_dim: int, hidden: int = 32) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Linear(hidden, 4),
        )

    def forward(self, features: torch.Tensor) -> torch.Tensor:
        return self.net(features)


@dataclass
class Surrogate:
    """Standardized cost network. ``predict`` returns decoded measurements."""

    net: CostNet
    x_mean: np.ndarray
    x_std: np.ndarray
    y_mean: np.ndarray
    y_std: np.ndarray
    history: dict[str, float]

    def predict_raw(self, features: np.ndarray) -> np.ndarray:
        """Network outputs in the encoded target space (log bytes, and so on)."""

        self.net.eval()
        x = (features - self.x_mean) / self.x_std
        with torch.no_grad():
            out = self.net(torch.tensor(x, dtype=torch.float32)).numpy()
        return out * self.y_std + self.y_mean

    def predict(self, features: np.ndarray) -> np.ndarray:
        """Columns: bytes, remote fraction, shuffle bytes, seconds."""

        if features.size == 0:
            return np.zeros((0, 4))
        return decode_targets(self.predict_raw(features))


def _standardize(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = matrix.mean(axis=0)
    std = matrix.std(axis=0)
    std = np.where(std < 1e-6, 1.0, std)
    return (matrix - mean) / std, mean, std


def train_surrogate(
    features: np.ndarray,
    targets: np.ndarray,
    hidden: int = 32,
    epochs: int = 80,
    seed: int = 0,
    lr: float = 1e-3,
    patience: int = 12,
    val_fraction: float = 0.2,
) -> Surrogate:
    """Fit the cost network with early stopping on a row holdout.

    ``targets`` is the encoded four-column matrix from ``stack_targets``.
    """

    if features.ndim != 2 or targets.shape[0] != features.shape[0]:
        raise ValueError("features and targets are not aligned")
    if features.shape[0] < 4:
        raise ValueError("need at least 4 measured rows to fit the cost network")
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    x_std, x_mean, x_scale = _standardize(features.astype(np.float64))
    y_std, y_mean, y_scale = _standardize(targets.astype(np.float64))
    n = x_std.shape[0]
    order = rng.permutation(n)
    n_val = int(round(n * val_fraction)) if n >= 10 else max(1, n // 5)
    n_val = min(n_val, n - 3)
    val_idx = order[:n_val]
    fit_idx = order[n_val:]
    x_fit = torch.tensor(x_std[fit_idx], dtype=torch.float32)
    y_fit = torch.tensor(y_std[fit_idx], dtype=torch.float32)
    x_val = torch.tensor(x_std[val_idx], dtype=torch.float32)
    y_val = torch.tensor(y_std[val_idx], dtype=torch.float32)
    net = CostNet(features.shape[1], hidden)
    optimizer = torch.optim.Adam(net.parameters(), lr=lr, weight_decay=1e-3)
    best_state = {k: v.detach().clone() for k, v in net.state_dict().items()}
    best_val = float("inf")
    stale = 0
    for _epoch in range(epochs):
        net.train()
        optimizer.zero_grad()
        loss = torch.mean((net(x_fit) - y_fit) ** 2)
        loss.backward()
        optimizer.step()
        net.eval()
        with torch.no_grad():
            val = float(torch.mean((net(x_val) - y_val) ** 2))
        if val + 1e-6 < best_val:
            best_val = val
            best_state = {k: v.detach().clone() for k, v in net.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= patience:
                break
    net.load_state_dict(best_state)
    net.eval()
    return Surrogate(
        net=net,
        x_mean=x_mean,
        x_std=x_scale,
        y_mean=y_mean,
        y_std=y_scale,
        history={"val_mse": best_val, "rows": float(n), "val_rows": float(n_val)},
    )
