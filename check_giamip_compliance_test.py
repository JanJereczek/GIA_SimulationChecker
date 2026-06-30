#!/usr/bin/env python3
"""Run the GIAMIP compliance checker on output/CUBoulder-SemiAnalytic."""

import sys
import os

_CONDA_PYTHON = os.path.expanduser("~/.miniconda3/envs/isschecker/bin/python")
if os.path.exists(_CONDA_PYTHON) and os.path.realpath(sys.executable) != os.path.realpath(_CONDA_PYTHON):
    os.execv(_CONDA_PYTHON, [_CONDA_PYTHON] + sys.argv)

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import giamip_compliance_checker as checker

SOURCE_PATH = REPO_ROOT / "output" / "CUBoulder-SemiAnalytic"
FORCING_PATH = REPO_ROOT / "input" / "iceHistory-PaleoMIST_1a.nc"

summary = checker.run_checker(
    source_path=str(SOURCE_PATH),
    workdir=str(REPO_ROOT),
    forcing_path=str(FORCING_PATH),
)

print(f"\nTotal errors                 : {summary['total_errors']}")
print(f"  Missing mandatory variable : {summary['total_file_errors']}")
print(f"  Naming                     : {summary['total_naming_errors']}")
print(f"  Numerical                  : {summary['total_num_errors']}")
print(f"  Spatial                    : {summary['total_spatial_errors']}")
print(f"  Dimension                  : {summary['total_dim_errors']}")
print(f"  Time                       : {summary['total_time_errors']}")
print(f"  Metadata attributes        : {summary['total_meta_attr_errors']}")
print(f"  Coordinate attributes      : {summary['total_coord_attr_errors']}")
print(f"  Variable attributes        : {summary['total_var_attr_errors']}")
print(f"\nFull log                     : {summary['log_path']}")
