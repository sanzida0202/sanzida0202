"""Command-line entry point for the MSI vs. IHC concordance analysis."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from msi_ihc_concordance.analysis import run_analysis


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Microsatellite instability vs. immunohistochemistry concordance analysis",
    )
    parser.add_argument(
        "--msi-data",
        type=Path,
        default=Path("data/msi_calls.csv"),
        help="Path to the MSI calls CSV file.",
    )
    parser.add_argument(
        "--ihc-data",
        type=Path,
        default=Path("data/ihc_results.csv"),
        help="Path to the IHC results CSV file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="Directory where analysis outputs will be written.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    summary = run_analysis(args.msi_data, args.ihc_data, args.output_dir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
