#!/usr/bin/env python3
"""
NetCDF file generator for GIAMIP GIA simulation output.

Generates GIAMIP-compliant NetCDF files with synthetic data on the prescribed
global Gaussian grid (257 lat × 513 lon nodes).  One file is written per variable,
following the naming convention:

    {variable_name}_{experiment_id}_{group_name}_{model_name}.nc

Time encoding:  "days since 0001-01-01 00:00:00", calendar "proleptic_gregorian".
All data variables use float32.
"""

import argparse
from pathlib import Path

import cftime
import netCDF4
import numpy as np
import xarray as xr

from giamip_compliance_checker import GIAMIP_VARIABLES, GIAMIP_VAR_META

# Full GIAMIP output grid
NLAT = 257
NLON = 513
TIME_UNITS = "days since 0001-01-01 00:00:00"
CALENDAR = "proleptic_gregorian"


def _make_grid():
    # Gaussian grids never include the exact poles; shift by half a cell.
    step = 180.0 / NLAT
    lat = np.linspace(-90.0 + step / 2, 90.0 - step / 2, NLAT, dtype=np.float32)
    lon = np.linspace(0.0, 360.0, NLON, endpoint=False, dtype=np.float32)
    return lat, lon


def _year_to_days(year: int) -> float:
    """Convert an integer year to days since 0001-01-01 (proleptic_gregorian)."""
    date = cftime.DatetimeProlepticGregorian(year, 1, 1)
    return float(cftime.date2num(date, TIME_UNITS, CALENDAR))


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

    # --- Build time axis ---
    year_values = np.linspace(start_year, end_year, n_steps, dtype=float)
    time_days = np.array([_year_to_days(int(y)) for y in year_values], dtype=np.float32)

    # --- Build data ---
    if dim_type == "lat_lon_t":
        lat, lon = _make_grid()
        data = np.random.uniform(0.1, 0.9, (n_steps, NLAT, NLON)).astype(np.float32)
        # For mask variables ensure [0, 1]
        if variable_name not in ("ocean_area_fraction", "land_ice_area_fraction"):
            data = data * 100.0  # give some spread for non-mask vars

        coords = {
            "time": ("time", time_days, {
                "units": TIME_UNITS, "calendar": CALENDAR,
            }),
            "lat": ("lat", lat, {"units": "degrees_north", "long_name": "latitude"}),
            "lon": ("lon", lon, {"units": "degrees_east", "long_name": "longitude"}),
        }
        data_vars = {
            variable_name: (
                ("time", "lat", "lon"),
                data,
                _var_attrs(var_meta),
            )
        }

    elif dim_type == "t":
        data = np.random.uniform(1e10, 1e12, n_steps).astype(np.float32)

        coords = {
            "time": ("time", time_days, {
                "units": TIME_UNITS, "calendar": CALENDAR,
            }),
        }
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
        data = np.zeros((n, n), dtype=np.float32)
        # Fill upper-triangular entries (valid degree/order combinations)
        for d in range(n):
            data[d, : d + 1] = np.random.uniform(-1e-6, 1e-6, d + 1).astype(np.float32)

        coords = {
            "degree": ("degree", np.arange(n, dtype=np.int32)),
            "order": ("order", np.arange(n, dtype=np.int32)),
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
        "contact_name": "Synthetic Test",
        "contact_email": "test@example.org",
        "reference_frame": "CM",
    })

    filename = f"{variable_name}_{experiment_id}_{group}_{model}.nc"
    output_path = output_dir / filename

    unlimited = ("time",) if dim_type in ("lat_lon_t", "t") else ()
    encoding = {}
    if variable_name in ds.data_vars:
        encoding[variable_name] = {"dtype": "f4", "_FillValue": None}
    if "time" in ds.coords:
        encoding["time"] = {"dtype": "f4", "_FillValue": None}

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
        default="./models/GIAMIP/Exp01/CORE",
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
