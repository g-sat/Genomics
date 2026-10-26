"""Figure for the chromosome 22 rule check. Reads results/real_slice/measurements.json."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
payload = json.loads((ROOT / "results" / "real_slice" / "measurements.json").read_text(encoding="utf-8"))
reads = payload["mean_read_bytes"]
fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))
width = 0.36
selective_labels = ["Inside one family", "Family plus 12 others"]
selective_family = [reads["family_inside"], reads["family_cross"]]
selective_all = [reads["all_samples_window"], reads["all_samples_window"]]
x = np.arange(2)
axes[0].bar(x - width / 2, np.array(selective_all) / 1024, width, color="#4c78a8", label="All-samples tile")
axes[0].bar(x + width / 2, np.array(selective_family) / 1024, width, color="#f58518", label="Family files")
axes[0].set_xticks(x)
axes[0].set_xticklabels(selective_labels)
axes[0].set_ylabel("Opaque bytes read (KB)")
axes[0].set_title("A check over 200 sites")
axes[0].legend(frameon=False, fontsize=8)
scan_labels = ["All-samples", "Family files", "Populations"]
scan_vals = [reads["scan_all"], reads["scan_family"], reads["scan_population"]]
axes[1].bar(np.arange(3), np.array(scan_vals) / 1e6, color=["#4c78a8", "#f58518", "#54a24b"])
axes[1].set_xticks(np.arange(3))
axes[1].set_xticklabels(scan_labels)
axes[1].set_ylabel("Opaque bytes read (MB)")
axes[1].set_title("A scan of all 3,000 sites")
fig.suptitle("Chromosome 22, Phase 3 dosages, 64 KB blocks", y=1.02)
fig.tight_layout()
path = ROOT / "results" / "real_slice" / "rule.png"
fig.savefig(path, dpi=140)
(ROOT / "paper" / "figures" / "real_slice.png").write_bytes(path.read_bytes())
print(path)
