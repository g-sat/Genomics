"""Download the 1000 Genomes phase 3 pedigree and keep complete parent pairs.

Dosages are still simulated later. People with one recorded parent are dropped,
because the genotype simulator requires both parents or neither. A child is
dropped again if a parent was removed. The written file is PLINK PED:
family, individual, father, mother, sex, phenotype.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped"
RAW = ROOT / "data" / "pedigrees" / "20130606_g1k.ped"
CLEAN = ROOT / "data" / "pedigrees" / "1kg_complete_parents.ped"
MISSING = {"0", "-9", ""}


def _rows(text: str) -> list[tuple[str, str, str, str, str]]:
    parsed: list[tuple[str, str, str, str, str]] = []
    for line in text.splitlines():
        raw = line.strip()
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split()
        if len(parts) < 5:
            continue
        fam, iid, father, mother, sex = parts[:5]
        if sex not in {"1", "2"}:
            continue
        parsed.append((fam, iid, father, mother, sex))
    return parsed


def _complete(rows: list[tuple[str, str, str, str, str]]) -> list[tuple[str, str, str, str, str]]:
    current = list(rows)
    changed = True
    while changed:
        changed = False
        known = {iid for _fam, iid, _father, _mother, _sex in current}
        kept = []
        for fam, iid, father, mother, sex in current:
            father_known = father not in MISSING and father in known and father != iid
            mother_known = mother not in MISSING and mother in known and mother != iid
            if father_known == mother_known:
                if father_known:
                    kept.append((fam, iid, father, mother, sex))
                else:
                    kept.append((fam, iid, "0", "0", sex))
            else:
                changed = True
        if len(kept) != len(current):
            changed = True
        current = kept
    return current


def main() -> None:
    RAW.parent.mkdir(parents=True, exist_ok=True)
    if not RAW.exists():
        print(f"downloading {SOURCE}", flush=True)
        urllib.request.urlretrieve(SOURCE, RAW)
    rows = _rows(RAW.read_text(encoding="utf-8"))
    if len({iid for _fam, iid, _f, _m, _s in rows}) != len(rows):
        raise RuntimeError("1000 Genomes PED has duplicate individual ids")
    kept = _complete(rows)
    trios = sum(1 for _fam, _iid, father, _mother, _sex in kept if father != "0")
    lines = [f"{fam} {iid} {father} {mother} {sex} -9" for fam, iid, father, mother, sex in kept]
    CLEAN.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"source rows {len(rows)}")
    print(f"kept {len(kept)}")
    print(f"non-founders {trios}")
    print(f"wrote {CLEAN}")


if __name__ == "__main__":
    main()
