"""Figure for where the pedigree gap holds and where it closes.

    python experiments/plot_boundary.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

ORDER = (
    ("select_05", "5% of families"),
    ("select_15", "15% of families"),
    ("mix_50", "50/50 scan+family"),
    ("select_40", "40% of families"),
    ("scale", "102–147 samples"),
    ("overlap", "3 generations"),
)
SEEDS = (7, 11, 19)


def _lookup(summary: list[dict]) -> dict[str, dict[int, dict]]:
    out: dict[str, dict[int, dict]] = {}
    for block in summary:
        out[block["name"]] = {int(row["seed"]): row for row in block["rows"]}
    return out


def main() -> None:
    summary = json.loads((ROOT / "results" / "stress" / "summary.json").read_text(encoding="utf-8"))
    rows = _lookup(summary)
    labels = [label for _name, label in ORDER]
    objective = {seed: [] for seed in SEEDS}
    heldout = {seed: [] for seed in SEEDS}
    for name, _label in ORDER:
        for seed in SEEDS:
            row = rows[name][seed]
            objective[seed].append(row["blind_over_pedpack"])
            heldout[seed].append(row["heldout_bytes_blind"] / row["heldout_bytes_pedpack"])

    x = np.arange(len(labels))
    width = 0.24
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4))

    ax = axes[0]
    for i, seed in enumerate(SEEDS):
        ax.bar(x + (i - 1) * width, objective[seed], width, label=f"seed {seed}")
    ax.axhline(1.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Fit objective, no-pedigree / PedPack")
    ax.set_title("When the pedigree ranker is cheaper")
    ax.set_ylim(0.95, 1.28)
    ax.legend(frameon=False, ncol=3)

    ax = axes[1]
    for i, seed in enumerate(SEEDS):
        ax.bar(x + (i - 1) * width, heldout[seed], width, label=f"seed {seed}")
    ax.axhline(1.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Held-out family bytes, no-pedigree / PedPack")
    ax.set_title("Bytes on families the search did not fit")
    ax.legend(frameon=False, ncol=3)

    fig.tight_layout()
    path = ROOT / "results" / "stress" / "boundary.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print(path)


if __name__ == "__main__":
    main()
