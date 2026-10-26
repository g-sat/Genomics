"""Figure for the four synthetic layout extensions."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BLUE = "#4c78a8"
ORANGE = "#f58518"
GREEN = "#54a24b"
RED = "#e45756"

fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.4))

ax = axes[0, 0]
labels = ["One-family\ncopy", "Sixteen split\nfamilies", "Full\nrebuild", "Log after\ncompaction"]
values = [120, 1920, 4800, 8000]
colors = [ORANGE, ORANGE, BLUE, GREEN]
ax.bar(np.arange(4), values, color=colors)
ax.set_xticks(np.arange(4))
ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("Bytes")
ax.set_title("A. Append-only log, 120 people")
ax.set_ylim(0, 9200)

ax = axes[0, 1]
scenarios = ["Trio in\none block", "Two trios,\nshared parents", "Trio across\nthree blocks"]
block = np.array([65536, 65536, 196608])
stitch = np.array([600, 800, 600])
x = np.arange(3)
width = 0.36
ax.bar(x - width / 2, block, width, color=BLUE, label="Rounded blocks")
ax.bar(x + width / 2, stitch, width, color=ORANGE, label="Pointer stitch")
ax.set_xticks(x)
ax.set_xticklabels(scenarios, fontsize=8)
ax.set_ylabel("Bytes")
ax.set_yscale("log")
ax.set_title("B. Pointer stitch, 64 KB blocks")
ax.legend(frameon=False, fontsize=8)

ax = axes[1, 0]
h = np.linspace(0, 1, 201)
ped = 1080 - 720 * h
phe = 360 + 720 * h
split = 960 - 480 * h
ax.plot(h, ped, color=BLUE, label="Pedigree tile")
ax.plot(h, phe, color=ORANGE, label="Phenotype tile")
ax.plot(h, split, color=GREEN, label="Split files")
ax.axvline(0.5, color="#888888", linewidth=0.8, linestyle="--")
ax.set_xlabel("Mendelian fraction")
ax.set_ylabel("Bytes per joint query")
ax.set_title("C. Joint genotype and phenotype")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1200)
ax.legend(frameon=False, fontsize=8)

ax = axes[1, 1]
acc_labels = ["Proof", "One trio\nreturned", "All 32\npayloads"]
acc_values = [192, 120, 3840]
ax.bar(np.arange(3), acc_values, color=[GREEN, ORANGE, BLUE])
ax.set_xticks(np.arange(3))
ax.set_xticklabels(acc_labels, fontsize=8)
ax.set_ylabel("Bytes")
ax.set_title("D. Trio accumulator")
ax.set_ylim(0, 4600)

fig.tight_layout()
out = ROOT / "results" / "extensions" / "extensions.png"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=160)
dest = ROOT / "paper" / "figures" / "extensions.png"
dest.write_bytes(out.read_bytes())
print(out)
