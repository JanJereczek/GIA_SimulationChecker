#!/usr/bin/env python3
"""
NetCDF file generator for GIAMIP GIA simulation output.

Generates GIAMIP-compliant NetCDF files with synthetic data on the prescribed
global Gaussian grid (257 lat × 513 lon nodes).  One file is written per variable,
following the naming convention:

    {variable_name}_{experiment_id}_{group_name}_{model_name}.nc

Time encoding:  'time' holds days since 2000-01-01 on a 365-day calendar, and the
companion 'year' variable holds the year of each step (time/365).
Each data variable uses the precision requested by the data request (float32 for
the gridded fields, float64 for the scalar time series and Clm/Slm).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr

# Allow running this script directly (python generate/generate_giamip_test_files.py):
# ensure the repository root (which holds giamip_compliance_checker.py) is importable.
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from giamip_compliance_checker import (
    GIAMIP_COORD_LONG_NAMES,
    GIAMIP_DAYS_PER_YEAR,
    GIAMIP_TIME_CALENDAR,
    GIAMIP_TIME_UNITS,
    GIAMIP_VARIABLES,
    GIAMIP_VAR_META,
    GIAMIP_YEAR_VAR,
)

# Full GIAMIP output grid
NLAT = 257
NLON = 513
TIME_UNITS = GIAMIP_TIME_UNITS

# Coordinate attributes required by the data request.
COORD_ATTRS = {
    "time": {"units": TIME_UNITS, "calendar": GIAMIP_TIME_CALENDAR,
             "long_name": GIAMIP_COORD_LONG_NAMES["time"], "standard_name": "time"},
    "year": {"units": "year", "long_name": GIAMIP_COORD_LONG_NAMES["year"]},
    "lat": {"units": "degrees_north", "long_name": "Latitude",
            "standard_name": "latitude"},
    "lon": {"units": "degrees_east", "long_name": "Longitude",
            "standard_name": "longitude"},
    "degree": {"units": "1", "long_name": "Spherical harmonic degree"},
    "order": {"units": "1", "long_name": "Spherical harmonic order"},
}


def _time_coords(years: np.ndarray) -> dict:
    """Build the 'time' (days since 2000-01-01) and 'year' coordinate pair."""
    years = np.asarray(years, dtype=np.int32)
    days = years.astype(np.float64) * GIAMIP_DAYS_PER_YEAR
    return {
        "time": ("time", days, COORD_ATTRS["time"]),
        GIAMIP_YEAR_VAR: ("time", years, COORD_ATTRS["year"]),
    }


def _make_grid():
    # Gaussian grids never include the exact poles; shift by half a cell.
    step = 180.0 / NLAT
    lat = np.linspace(-90.0 + step / 2, 90.0 - step / 2, NLAT, dtype=np.float64)
    lon = np.linspace(0.0, 360.0, NLON, endpoint=False, dtype=np.float64)
    return lat, lon


def create_giamip_file(
    output_dir: Path,
    variable_name: str,
    experiment_id: str,
    group: str,
    model: str,
    start_year: int = 1,
    end_year: int = 1001,
    n_steps: int = 2,
) -> Path:
    """
    Write one GIAMIP-compliant NetCDF file for *variable_name*.

    Returns the Path of the created file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    var_meta = GIAMIP_VAR_META[variable_name]
    dim_type = var_meta["dim"]   # "lat_lon_t" | "t" | "degree_order"
    precision = var_meta["precision"]   # np.float32 | np.float64

    # --- Build time axis (days since 2000-01-01, plus the year of each step) ---
    time_years = np.linspace(start_year, end_year, n_steps).round().astype(np.int32)

    # --- Build data ---
    if dim_type == "lat_lon_t":
        lat, lon = _make_grid()
        data = np.random.uniform(0.1, 0.9, (n_steps, NLAT, NLON)).astype(precision)
        # For mask variables ensure [0, 1]
        if variable_name not in ("ocean_area_fraction", "land_ice_area_fraction"):
            data = data * 100.0  # give some spread for non-mask vars

        coords = {
            **_time_coords(time_years),
            "lat": ("lat", lat, COORD_ATTRS["lat"]),
            "lon": ("lon", lon, COORD_ATTRS["lon"]),
        }
        data_vars = {
            variable_name: (
                ("time", "lat", "lon"),
                data,
                _var_attrs(var_meta),
            )
        }

    elif dim_type == "t":
        # Draw values comfortably inside the variable's plausible bounds.
        lo, hi = var_meta.get("bounds", (0.0, 1.0))
        mid, half = 0.5 * (lo + hi), 0.2 * (hi - lo)
        data = np.random.uniform(mid - half, mid + half, n_steps).astype(precision)

        coords = dict(_time_coords(time_years))
        data_vars = {
            variable_name: (
                ("time",),
                data,
                _var_attrs(var_meta),
            )
        }

    elif dim_type == "degree_order":
        # Stokes coefficients up to degree 96 stored as (97, 97) matrix
        max_degree = 96
        n = max_degree + 1
        data = np.zeros((n, n), dtype=precision)
        # Fill upper-triangular entries (valid degree/order combinations)
        for d in range(n):
            data[d, : d + 1] = np.random.uniform(-1e-6, 1e-6, d + 1).astype(precision)

        coords = {
            "degree": ("degree", np.arange(n, dtype=np.float64), COORD_ATTRS["degree"]),
            "order": ("order", np.arange(n, dtype=np.float64), COORD_ATTRS["order"]),
        }
        data_vars = {
            variable_name: (
                ("degree", "order"),
                data,
                _var_attrs(var_meta),
            )
        }

    else:
        raise ValueError(f"Unknown dim_type '{dim_type}' for variable '{variable_name}'")

    ds = xr.Dataset(data_vars, coords=coords)
    ds.attrs.update({
        "title": f"GIAMIP synthetic data - {variable_name}",
        "Conventions": "CF-1.7",
        "group": group,
        "model": model,
        "contributors_name": "Synthetic Test",
        "contributors_email": "test@example.org",
    })

    filename = f"{variable_name}_{experiment_id}_{group}_{model}.nc"
    output_path = output_dir / filename

    unlimited = ("time",) if dim_type in ("lat_lon_t", "t") else ()
    encoding = {}
    if variable_name in ds.data_vars:
        encoding[variable_name] = {"dtype": np.dtype(precision).str[1:],
                                   "_FillValue": None}
    if "time" in ds.coords:
        encoding["time"] = {"dtype": "f8", "_FillValue": None}
    if GIAMIP_YEAR_VAR in ds.coords:
        encoding[GIAMIP_YEAR_VAR] = {"dtype": "i4", "_FillValue": None}

    ds.to_netcdf(output_path, unlimited_dims=unlimited, encoding=encoding)
    return output_path


def _var_attrs(meta: dict) -> dict:
    attrs = {
        "long_name": meta["long_name"],
        "units": meta["units"],
    }
    if meta.get("standard_name") is not None:
        attrs["standard_name"] = meta["standard_name"]
    return attrs


def create_all_mandatory_files(
    output_dir: Path,
    experiment_id: str = "Exp01",
    group: str = "TESTGROUP",
    model: str = "TESTMODEL",
    start_year: int = 1,
    end_year: int = 1001,
    n_steps: int = 2,
) -> list[Path]:
    """Create one file for every mandatory GIAMIP variable."""
    mandatory = [v["variable"] for v in GIAMIP_VARIABLES if v["mandatory"]]
    created = []
    for var in mandatory:
        path = create_giamip_file(
            output_dir=output_dir,
            variable_name=var,
            experiment_id=experiment_id,
            group=group,
            model=model,
            start_year=start_year,
            end_year=end_year,
            n_steps=n_steps,
        )
        created.append(path)
    print(f"Created {len(created)} mandatory GIAMIP files in {output_dir}")
    return created


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate GIAMIP-compliant synthetic NetCDF files."
    )
    parser.add_argument(
        "--output-dir",
        default="./output/GIAMIP/Exp01/CORE",
        help="Directory for output files.",
    )
    parser.add_argument(
        "--experiment-id", default="Exp01",
        help="Experiment identifier (default: Exp01).",
    )
    parser.add_argument(
        "--group", default="TESTGROUP",
        help="Group name (default: TESTGROUP).",
    )
    parser.add_argument(
        "--model", default="TESTMODEL",
        help="Model name (default: TESTMODEL).",
    )
    parser.add_argument(
        "--start-year", type=int, default=1,
        help="First year of the time axis (default: 1).",
    )
    parser.add_argument(
        "--end-year", type=int, default=1001,
        help="Last year of the time axis (default: 1001).",
    )
    parser.add_argument(
        "--n-steps", type=int, default=2,
        help="Number of time steps (default: 2).",
    )
    parser.add_argument(
        "--variable",
        choices=[v["variable"] for v in GIAMIP_VARIABLES],
        help="Generate a single variable only (default: all mandatory).",
    )
    args = parser.parse_args()

    if args.variable:
        create_giamip_file(
            output_dir=Path(args.output_dir),
            variable_name=args.variable,
            experiment_id=args.experiment_id,
            group=args.group,
            model=args.model,
            start_year=args.start_year,
            end_year=args.end_year,
            n_steps=args.n_steps,
        )
    else:
        create_all_mandatory_files(
            output_dir=Path(args.output_dir),
            experiment_id=args.experiment_id,
            group=args.group,
            model=args.model,
            start_year=args.start_year,
            end_year=args.end_year,
            n_steps=args.n_steps,
        )
