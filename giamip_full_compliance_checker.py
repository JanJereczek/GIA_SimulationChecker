#!/usr/bin/env python3
"""Check all submitted GIAMIP models at once.

For every model directory under --source-dir, loops over its experiment
subdirectories and runs the per-experiment compliance checker, selecting each
experiment's forcing file from --forcing-dir (see experiments.md).

Example:
    python giamip_full_compliance_checker.py \\
        --source-dir ./output \\
        --forcing-dir ./input
"""

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import giamip_compliance_checker as checker


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check all submitted GIAMIP models at once."
    )
    parser.add_argument(
        "--source-dir",
        required=True,
        help="Top-level output directory containing model subdirectories"
             " (e.g. ./output).",
    )
    parser.add_argument(
        "--forcing-dir",
        required=True,
        help="Directory containing the forcing NetCDF files (e.g. ./input).",
    )
    args = parser.parse_args()

    checker.warn_if_wrong_environment()
    checker.run_full_checker(
        output_dir=args.source_dir,
        forcing_dir=args.forcing_dir,
    )


if __name__ == "__main__":
    main()
