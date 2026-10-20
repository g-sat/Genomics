"""Run a PedPack simulator study.

Example:
    python experiments/run_study.py --config configs/small.json --out results/small
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pedpack.report import write_report  # noqa: E402
from pedpack.study import run_study  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare pedigree-aware layouts on the block simulator")
    parser.add_argument("--config", type=Path, default=ROOT / "configs" / "small.json")
    parser.add_argument("--out", type=Path, default=ROOT / "results" / "small")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    report = run_study(config, verbose=not args.quiet)
    write_report(report, args.out)
    print(f"wrote {args.out / 'report.md'}")


if __name__ == "__main__":
    main()
