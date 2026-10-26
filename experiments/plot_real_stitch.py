"""Figure for the chromosome 22 pointer stitch and Parquet row groups."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
record = json.loads((ROOT / "results" / "real_stitch" / "measurements.json").read_text(encoding="utf-8"))
full = record["sites_3000"]
parquet = record["parquet"]
BLUE = "#4c78a8"
ORANGE = "#f58518"
GREEN = "#54a24b"

fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))
labels = ["Stitch", "One block", "Packed file"]
bytes_ = [full["stitch_bytes"], full["rounded_bytes"], full["cohort_file_bytes"]]
seconds = [
    full["unbuffered_stitch_seconds"] * 1e3,
    full["unbuffered_trio_seconds"] * 1e3,
    full["unbuffered_cohort_seconds"] * 1e3,
]
x = np.arange(3)
axes[0].bar(x, bytes_, color=[ORANGE, BLUE, GREEN])
axes[0].set_xticks(x)
axes[0].set_xticklabels(labels)
axes[0].set_yscale("log")
axes[0].set_ylabel("Bytes")
axes[0].set_title("Unbuffered read, 3,000 sites")
axes[1].bar(
    np.arange(2),
    [parquet["trio_row_group_bytes"], parquet["whole_file_bytes"]],
    color=[ORANGE, BLUE],
)
axes[1].set_xticks(np.arange(2))
axes[1].set_xticklabels(["Trio row group", "All-samples file"])
axes[1].set_yscale("log")
axes[1].set_ylabel("Parquet bytes")
axes[1].set_title("Same dosages, one column of bytes")
fig.tight_layout()
out = ROOT / "results" / "real_stitch" / "stitch.png"
fig.savefig(out, dpi=160)
dest = ROOT / "paper" / "figures" / "real_stitch.png"
dest.write_bytes(out.read_bytes())
print(out)
