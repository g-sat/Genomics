"""Zstandard size of nuclear-family tiles against unrelated tiles of the same shape.

The cost model shrinks a tile by ``compress_base * (1 - kinship_alpha * mean kinship)``.
This script measures a generic codec on the simulated haplotypes and dosages.
It is not a PBWT implementation. A gap here means related columns are smaller
under Zstandard. A tie means this codec does not see the kinship term.

    python experiments/measure_compression.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pedpack.cost import payload_bytes  # noqa: E402
from pedpack.genotype import build_cohort  # noqa: E402
from pedpack.grouping import nuclear_groups  # noqa: E402
from pedpack.layout import group_mean_kinship  # noqa: E402
from pedpack.pedigree import simulate_pedigree  # noqa: E402

WINDOW = 800
COMPRESS_LEVEL = 3
BYTES_PER_GENOTYPE = 400.0
COMPRESS_BASE = 0.45
KINSHIP_ALPHA = 0.55


def _zstd_size(raw: bytes) -> int:
    import pyarrow as pa

    return len(pa.Codec("zstd", compression_level=COMPRESS_LEVEL).compress(raw))


def _tile(matrix: np.ndarray, columns: list[int], v0: int, v1: int) -> bytes:
    block = np.ascontiguousarray(matrix[v0:v1][:, columns])
    return block.tobytes()


def _haplotype_tile(
    paternal: np.ndarray,
    maternal: np.ndarray,
    columns: list[int],
    v0: int,
    v1: int,
) -> bytes:
    pat = paternal[v0:v1][:, columns]
    mat = maternal[v0:v1][:, columns]
    stacked = np.stack([pat, mat], axis=-1)
    return np.ascontiguousarray(stacked).tobytes()


def _windows(n_variants: int) -> list[tuple[int, int]]:
    if n_variants < WINDOW:
        return [(0, n_variants)]
    return [(v0, v0 + WINDOW) for v0 in range(0, n_variants - WINDOW + 1, WINDOW)]


def _paired(columns: list[int], founders: np.ndarray, rng: np.random.Generator) -> list[int] | None:
    if int(founders.size) < len(columns):
        return None
    picked = rng.choice(founders, size=len(columns), replace=False)
    return [int(x) for x in picked.tolist()]


def _sum_sizes(builder, columns: list[int], windows: list[tuple[int, int]]) -> tuple[int, int]:
    raw = 0
    compressed = 0
    for v0, v1 in windows:
        blob = builder(columns, v0, v1)
        raw += len(blob)
        compressed += _zstd_size(blob)
    return raw, compressed


def _regime(name: str, config: dict, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    pedigree = simulate_pedigree(
        n_families=int(config["n_families"]),
        generations=int(config["generations"]),
        mean_children=float(config["mean_children"]),
        cross_rate=float(config["cross_rate"]),
        rng=rng,
        min_people=int(config.get("min_people", 24)),
    )
    cohort, paternal, maternal = build_cohort(
        pedigree,
        n_variants=int(config["n_variants"]),
        rng=rng,
        bp_spacing=int(config.get("bp_spacing", 100)),
        return_haplotypes=True,
    )
    founders = np.flatnonzero((cohort.father < 0) & (cohort.mother < 0))
    windows = _windows(cohort.n_variants)
    control_rng = np.random.default_rng(seed + 10007)
    rows = []
    for group in nuclear_groups(cohort):
        if len(group) < 3:
            continue
        control = _paired(group, founders, control_rng)
        if control is None:
            continue
        phi = group_mean_kinship(np.asarray(group, dtype=np.int32), cohort.kinship)
        phi_control = group_mean_kinship(np.asarray(control, dtype=np.int32), cohort.kinship)
        dosage_family = _sum_sizes(lambda cols, v0, v1: _tile(cohort.genotypes, cols, v0, v1), group, windows)
        dosage_control = _sum_sizes(lambda cols, v0, v1: _tile(cohort.genotypes, cols, v0, v1), control, windows)
        hap_family = _sum_sizes(
            lambda cols, v0, v1: _haplotype_tile(paternal, maternal, cols, v0, v1),
            group,
            windows,
        )
        hap_control = _sum_sizes(
            lambda cols, v0, v1: _haplotype_tile(paternal, maternal, cols, v0, v1),
            control,
            windows,
        )
        predicted_family = payload_bytes(
            len(group), WINDOW, BYTES_PER_GENOTYPE, phi, COMPRESS_BASE, KINSHIP_ALPHA
        )
        predicted_control = payload_bytes(
            len(group), WINDOW, BYTES_PER_GENOTYPE, 0.0, COMPRESS_BASE, KINSHIP_ALPHA
        )
        rows.append(
            {
                "group_size": len(group),
                "mean_kinship": phi,
                "control_mean_kinship": phi_control,
                "dosage_family_compressed": dosage_family[1],
                "dosage_control_compressed": dosage_control[1],
                "dosage_raw": dosage_family[0],
                "hap_family_compressed": hap_family[1],
                "hap_control_compressed": hap_control[1],
                "hap_raw": hap_family[0],
                "model_ratio": predicted_family / predicted_control,
            }
        )
    return {
        "regime": name,
        "seed": seed,
        "n_samples": cohort.n_samples,
        "n_variants": cohort.n_variants,
        "n_founders": int(founders.size),
        "n_groups": len(rows),
        "window": WINDOW,
        "groups": rows,
    }


def _ratio(num: np.ndarray, den: np.ndarray) -> float:
    total_den = float(den.sum())
    if total_den <= 0:
        return float("nan")
    return float(num.sum() / total_den)


def _summarize(seed_rows: list[dict]) -> dict:
    groups = [group for row in seed_rows for group in row["groups"]]
    if not groups:
        raise RuntimeError("no nuclear families large enough to compress")
    hap_family = np.asarray([g["hap_family_compressed"] for g in groups], dtype=np.float64)
    hap_control = np.asarray([g["hap_control_compressed"] for g in groups], dtype=np.float64)
    dos_family = np.asarray([g["dosage_family_compressed"] for g in groups], dtype=np.float64)
    dos_control = np.asarray([g["dosage_control_compressed"] for g in groups], dtype=np.float64)
    hap_raw = np.asarray([g["hap_raw"] for g in groups], dtype=np.float64)
    dos_raw = np.asarray([g["dosage_raw"] for g in groups], dtype=np.float64)
    model = np.asarray([g["model_ratio"] for g in groups], dtype=np.float64)
    kinship = np.asarray([g["mean_kinship"] for g in groups], dtype=np.float64)
    return {
        "n_groups": len(groups),
        "mean_kinship": float(kinship.mean()),
        "mean_group_size": float(np.mean([g["group_size"] for g in groups])),
        "haplotype_family_over_raw": _ratio(hap_family, hap_raw),
        "haplotype_control_over_raw": _ratio(hap_control, hap_raw),
        "haplotype_family_over_control": _ratio(hap_family, hap_control),
        "dosage_family_over_raw": _ratio(dos_family, dos_raw),
        "dosage_control_over_raw": _ratio(dos_control, dos_raw),
        "dosage_family_over_control": _ratio(dos_family, dos_control),
        "model_family_over_unrelated": float(np.mean(model)),
    }


def _fmt(value: float) -> str:
    return f"{value:.3f}"


def _report(payload: dict) -> str:
    lines = [
        "# Zstandard on nuclear-family tiles",
        "",
        "Simulated haplotypes and dosages, same generator and seeds as the studies.",
        f"Codec is Zstandard level {COMPRESS_LEVEL}. Each tile is {WINDOW} variants,",
        "variant-major, and the control is an equal-sized draw of founders.",
        "The model ratio is `payload(family kinship) / payload(kinship 0)` with",
        f"`compress_base={COMPRESS_BASE}` and `kinship_alpha={KINSHIP_ALPHA}`.",
        "This is not PBWT. A generic codec can miss a haplotype transform that",
        "PBWT would catch, and the other way around.",
        "",
        "| Regime | People | Groups | Mean kinship | Hap family/raw | Hap founder/raw | Hap family/founder | Dosage family/founder | Model family/unrelated |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for regime in payload["regimes"]:
        summary = regime["summary"]
        people = ", ".join(str(row["n_samples"]) for row in regime["seeds"])
        lines.append(
            "| {name} | {people} | {groups} | {kin} | {hf} | {hc} | {hr} | {dr} | {model} |".format(
                name=regime["name"],
                people=people,
                groups=summary["n_groups"],
                kin=_fmt(summary["mean_kinship"]),
                hf=_fmt(summary["haplotype_family_over_raw"]),
                hc=_fmt(summary["haplotype_control_over_raw"]),
                hr=_fmt(summary["haplotype_family_over_control"]),
                dr=_fmt(summary["dosage_family_over_control"]),
                model=_fmt(summary["model_family_over_unrelated"]),
            )
        )
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    out = ROOT / "results" / "compression"
    out.mkdir(parents=True, exist_ok=True)
    regimes = []
    for name, path in (("locked", "pedigree_gap.json"), ("scale", "scale.json")):
        config = json.loads((ROOT / "configs" / path).read_text(encoding="utf-8"))
        seed_rows = []
        for seed in config["seeds"]:
            print(f"{name} seed {seed}", flush=True)
            seed_rows.append(_regime(name, config, int(seed)))
        regimes.append({"name": name, "seeds": seed_rows, "summary": _summarize(seed_rows)})
    payload = {
        "codec": "zstd",
        "compression_level": COMPRESS_LEVEL,
        "window": WINDOW,
        "compress_base": COMPRESS_BASE,
        "kinship_alpha": KINSHIP_ALPHA,
        "regimes": regimes,
    }
    (out / "measurements.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    text = _report(payload)
    (out / "report.md").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
