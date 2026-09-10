## Differences to GIAMIP file checker

`giamip_standardize_grids.py` (reproduced below) is a post-processing script written independently of this repository, to prepare submitted output for analysis. It is a **fixer**, not a reporter: it mutates a dataset into a canonical analysis form and writes `<var>_<model>_std.nc`, and anything it cannot fix raises and stops the run. Our checker, by contrast, opens files read-only and accumulates errors into a log. Its layout convention is `<modeldir>/<model>/<var>_<model>.nc` — no experiment id, no group name — so it operates downstream of submission, on files already sorted per model.

Note that the script is **not self-contained**: it imports `model_basename`, `model_dir`, `model_input_path` and `parse_models_arg` from a `plotting_utils` module that is not part of this repository, so it cannot be run from here as-is.

### What it does

The pipeline in `standardize_dataset`:

1. **`ensure_variable_name`** — if the canonical variable is absent, look for known aliases (`geoid_rate`, `geoid`, `dg`, `Geoid`, …) and case-insensitive matches, and rename.
2. **`standardize_dims_coords`** — rename time-coordinate aliases (`year`, `Year`, `Time`, `t`) to `time`, `latitude`/`longitude` dims and coords to `lat`/`lon`, transpose to `(lat, lon, time)`, then wrap negative longitudes into [0, 360] and re-sort.
3. **Precomputed-rate detection** — a variable with no `time` dimension (or carrying a `precomputed_rate` attribute) is treated as already being a rate; time handling and the time part of grid validation are skipped, and the attribute is stamped for downstream stages.
4. **`coerce_time_to_years`** — normalise the time coordinate to float64 years from int, `datetime64`, or cftime objects; for cftime with `has_year_zero`, subtract 2000.
5. **`validate_grid`** — 257 lat × 513 lon; longitude *regularly* spaced; latitude *irregularly* spaced (i.e. actually Gauss-Legendre).
6. **`convert_units`** — rescale values between m/dm/cm/mm.

`run_standardize_stage` is the programmatic entry point (called from a `run_pipeline.py`), alongside a CLI with `--var`, `--models`, `--modeldir`, `--input-unit`/`--target-unit`, `--suffix`, `--overwrite`, `--no-validate` and `--alias`.

### Checks it performs that we do not

Recorded for reference; the decision (2026-09-07) was that our checker already covers what it should, so none of these are currently implemented.

1. **Gaussian latitude spacing.** `validate_grid` fails if `lat` is *regularly* spaced. We check node counts, extents and monotonicity, and compare against the forcing grid when one is available — but with no forcing file we would pass an equally-spaced 257-point latitude axis. This also affects our own `generate/generate_giamip_test_files.py`, which builds latitude with `np.linspace`: the synthetic baseline is not on a Gaussian grid (the real PaleoMIST forcing has dlat varying between 0.6932° and 0.6990°) and would fail their validator.
2. **Longitude regular spacing.** Likewise never checked independently of the forcing file.
3. **Time axis must contain present day.** `infer_time_inc_present` requires an exact `time == 0.0` entry with a preceding step. We never check that year 0 (2000 CE) is on the axis, nor the data request's requirement that outputs include the simulation start and end timesteps — we get that only implicitly, and only when a forcing file is available.
4. **Non-year time encodings diagnosed properly.** They accept int, `datetime64` and cftime and convert them; we call `astype(int)` on whatever is there, so a `datetime64` axis becomes nanosecond integers and is reported as a forcing mismatch rather than "time is not stored as calendar years". Their `has_year_zero` → −2000 shift has no counterpart either. This is the same area as the known time-encoding divergence documented in [variables.md](variables.md).
5. **Unit *values*, not just the units string.** We compare the `units` attribute; they rescale m/dm/cm/mm. We have no check that magnitudes are consistent with the declared unit beyond the loose `bounds` (which would catch a mm-scale `delta_bed` only once it exceeded ±3000).
6. **Alias awareness in diagnostics.** They recognise `latitude`/`longitude`/`Time`/`t`/`year` coordinates and variable aliases such as `geoid`. Rejecting those is correct for a compliance checker, but we report "coordinate 'lat' not found" where we could say "found 'latitude' — rename to 'lat'".
7. **The precomputed-rate concept** (no time dimension). Analysis-specific, and out of scope for a compliance checker rather than a gap.

### Where their conventions conflict with the data request

- **Dimension order.** They transpose to `(lat, lon, time)`; the data request and our checker require `(time, lat, lon)`. Any `_std.nc` file fed back into our checker fails the dimension-ordering check.
- **`year` → `time`.** Their time-alias list renames a `year` coordinate to `time`. Version 3 of the data request defines `year` and `time` as two *distinct* required coordinate variables, so running this on a v3-conformant file would silently destroy that distinction.

Nothing in the script overlaps with our naming convention, mandatory-variable set, attribute checks, dtype/precision, fill values or bounds; that side is entirely ours.

## Code

#!/usr/bin/env python
"""Standardize GIAMIP model output NetCDFs to a common convention.

Callable as a Python module (see load_and_standardize / standardize_dataset)
or as a CLI that reads <var>_<model>.nc, standardizes it, and writes
<var>_<model><suffix>.nc (default suffix: _std).

Standardization steps applied to every dataset:
  1. Rename latitude/longitude/time aliases to lat/lon/time.
  2. Coerce time to float64 years (handles int, datetime64, cftime).
  3. Force lon into [0, 360] and transpose dims to (lat, lon, time).
  4. Optionally validate the grid is Gauss-Legendre 257x513
     (irregular lat, regular lon).
  5. Optionally convert the data-variable units.

Generic over variable name; the layout convention is:

    <modeldir>/<model>/<var>_<model>.nc            (raw input)
    <modeldir>/<model>/<var>_<model>_std.nc        (standardized output)

`<model>` may be passed either as a bare directory name (looked up under
`--modeldir`) or as a full path to a per-model directory.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Iterable

import numpy as np
import pandas as pd
import xarray as xr

from plotting_utils import (  # re-exported: parse_models_arg
    model_basename,
    model_dir,
    model_input_path,
    parse_models_arg,
)

__all__ = [
    "standardize_dims_coords",
    "coerce_time_to_years",
    "infer_time_inc_present",
    "validate_grid",
    "convert_units",
    "ensure_variable_name",
    "is_precomputed_rate",
    "standardize_dataset",
    "load_and_standardize",
    "parse_models_arg",
    "DEFAULT_VAR_ALIASES",
]

# xarray attr key set on the data variable when the standardizer determines the
# input file is already a rate of change (no time dimension to difference over).
PRECOMPUTED_RATE_ATTR = "precomputed_rate"


_UNIT_TO_METERS = {"m": 1.0, "dm": 10.0, "cm": 100.0, "mm": 1000.0}

# Known variable-name aliases per canonical GIAMIP variable. Some modelers name
# their fields differently on export; extend this dict when a new alias shows up
# in the wild, or override per-call via the `aliases` kwarg.
DEFAULT_VAR_ALIASES: dict[str, tuple[str, ...]] = {
    "delta_g": ("geoid_rate", "geoid", "delta_geoid", "dg", "Geoid", "GEOID"),
    "uplift":  ("uplift_rate", "vertical_displacement", "u_rate", "vlm"),
}


def ensure_variable_name(ds: xr.Dataset,
                         canonical: str,
                         aliases: tuple[str, ...] | list[str] | None = None
                         ) -> xr.Dataset:
    """Ensure `ds` has a data-variable named `canonical`.

    If it does, return `ds` unchanged. Otherwise search `aliases` (falls back
    to `DEFAULT_VAR_ALIASES[canonical]` when None) plus a case-insensitive
    match against all data variables; rename the first hit to `canonical`.

    Raises ValueError if none of the aliases (or a case-insensitive match) are
    present in the dataset.
    """
    if canonical in ds.data_vars:
        return ds

    candidates: list[str] = list(aliases) if aliases is not None else []
    candidates.extend(DEFAULT_VAR_ALIASES.get(canonical, ()))

    for name in candidates:
        if name in ds.data_vars:
            print(f"Renaming variable {name!r} -> {canonical!r}")
            return ds.rename({name: canonical})

    lower_map = {str(v).lower(): v for v in ds.data_vars}
    for name in [canonical, *candidates]:
        hit = lower_map.get(name.lower())
        if hit is not None and hit != canonical:
            print(f"Renaming variable {hit!r} -> {canonical!r} (case-insensitive match)")
            return ds.rename({hit: canonical})

    raise ValueError(
        f"Dataset has no variable named {canonical!r} or any known alias. "
        f"Available variables: {list(ds.data_vars)}. "
        f"Tried aliases: {candidates}. Extend DEFAULT_VAR_ALIASES or pass "
        f"aliases=... to standardize_dataset/load_and_standardize."
    )


_TIME_ALIASES = ('year', 'Year', 'YEAR', 'Time', 'TIME', 't')


def _rename_time_coord_alias(ds: xr.Dataset) -> xr.Dataset:
    """Rename a time-coord alias (year/Year/Time/t/...) to `time` if present.

    Idempotent: returns `ds` unchanged if `time` already exists or no alias is
    found. Callers that need `time` in the standardized position should invoke
    this before any dim-order operation.
    """
    if 'time' in ds.coords:
        return ds
    for alt in _TIME_ALIASES:
        if alt in ds.coords:
            print(f"Renaming coordinate {alt!r} to 'time'.")
            return ds.rename({alt: 'time'})
    return ds


def standardize_dims_coords(ds: xr.Dataset, data_var: str) -> xr.Dataset:
    """Rename latitude/longitude to lat/lon, force lon in [0,360], and put
    the target variable's dims into (lat, lon, time) order.

    Time-coord aliases (year, Time, t, ...) are renamed to `time` first so the
    transpose step sees the canonical dim name.
    """
    ds = _rename_time_coord_alias(ds)
    if 'latitude' in ds.dims:
        if 'lat' in ds.coords:
            for v in ds.data_vars:
                if 'latitude' in ds[v].dims:
                    ds = ds.swap_dims({'latitude': 'lat'})
        else:
            ds = ds.rename_dims({'latitude': 'lat'})

    if 'longitude' in ds.dims:
        if 'lon' in ds.coords:
            for v in ds.data_vars:
                if 'longitude' in ds[v].dims:
                    ds = ds.swap_dims({'longitude': 'lon'})
        else:
            ds = ds.rename_dims({'longitude': 'lon'})

    coord_renames = {}
    if 'latitude' in ds.coords and 'lat' not in ds.coords:
        coord_renames['latitude'] = 'lat'
    if 'longitude' in ds.coords and 'lon' not in ds.coords:
        coord_renames['longitude'] = 'lon'
    if coord_renames:
        ds = ds.rename(coord_renames)

    target_order = tuple(d for d in ('lat', 'lon', 'time') if d in ds[data_var].dims)
    if ds[data_var].dims != target_order:
        ds[data_var] = ds[data_var].transpose(*target_order)

    if 'lon' in ds.coords and np.any(ds['lon'].values < 0):
        lons = ds['lon'].values.copy()
        lons[lons < 0] = lons[lons < 0] + 360
        ds = ds.assign_coords(lon=lons).sortby('lon', ascending=True)
        print("Longitude values adjusted to be in the range [0, 360].\n")

    return ds


def coerce_time_to_years(ds: xr.Dataset) -> xr.Dataset:
    """Ensure `time` is a float64-years coord.

    Handles the same input types as calculate_EWH.py used to (int, datetime64,
    cftime with object dtype), plus common alias names (year/Year/Time/t).
    """
    ds = _rename_time_coord_alias(ds)
    if 'time' not in ds.coords:
        raise ValueError("Dataset missing required 'time' coordinate.")

    t = ds['time']
    if t.dtype == np.dtype('float64'):
        return ds

    if t.dtype in (np.dtype('int32'), np.dtype('int64')):
        print("Detected integer time coordinate; converting to float years.")
        return ds.assign_coords(time=t.values.astype('float64'))

    if np.issubdtype(t.dtype, np.datetime64):
        print("Detected numpy datetime64 time coordinate; converting to float years.")
        try:
            years = t.dt.year.values.astype('float64')
        except Exception:
            years = pd.to_datetime(t.values).year.to_numpy().astype('float64')
        return ds.assign_coords(time=years)

    if t.dtype == object:
        vals = t.values
        if len(vals) and hasattr(vals[0], 'year'):
            print("Detected cftime-like time coordinate; converting to float years.")
            years = np.array([getattr(v, 'year') for v in vals], dtype='float64')
            if getattr(vals[0], 'has_year_zero', False):
                print("Detected has_year_zero=True; shifting years by -2000.")
                years = years - 2000.0
            return ds.assign_coords(time=years)

    raise ValueError(f"Unsupported time dtype {t.dtype!r}; cannot coerce to float years.")


def infer_time_inc_present(ds: xr.Dataset) -> float:
    """Absolute time increment at present day (time == 0).

    Locates the index where time == 0.0 (within tolerance) and returns the
    absolute difference to the immediately preceding value.
    """
    if 'time' not in ds.coords:
        raise ValueError("Dataset missing 'time' coordinate; cannot infer increment.")
    t = ds['time'].values
    if t.size < 2:
        raise ValueError("Time coordinate must have at least 2 entries to infer increment.")

    idx0 = None
    for i, v in enumerate(t):
        try:
            if np.isclose(float(v), 0.0, atol=1e-9):
                idx0 = i
                break
        except Exception:
            pass
    if idx0 is None:
        raise ValueError("Could not locate time == 0.0 in time coordinate.")
    if idx0 == 0:
        raise ValueError("time == 0.0 is the first entry; no preceding timestep.")

    return abs(float(t[idx0]) - float(t[idx0 - 1]))


def is_precomputed_rate(ds: xr.Dataset, data_var: str) -> bool:
    """Return True when `ds[data_var]` looks like an already-computed rate.

    Currently: the variable has no `time` dimension, or the variable/dataset
    carries the `PRECOMPUTED_RATE_ATTR` attr set to a truthy value (1 or True;
    stored as int in NetCDF since netCDF4 doesn't accept bool attrs).
    """
    for attrs in (ds.attrs, ds[data_var].attrs if data_var in ds.data_vars else {}):
        val = attrs.get(PRECOMPUTED_RATE_ATTR)
        if val is None:
            continue
        # 1 / True / "true" all count. netCDF4 stores as int.
        if isinstance(val, str):
            if val.lower() in ('1', 'true', 'yes', 'on'):
                return True
        elif int(val) == 1:
            return True

    if data_var in ds.data_vars and 'time' not in ds[data_var].dims:
        return True
    return False


def validate_grid(ds: xr.Dataset,
                  expected_lat_len: int = 257,
                  expected_lon_len: int = 513,
                  *,
                  require_time: bool = True) -> None:
    """Assert the Gauss-Legendre 257x513 grid (irregular lat, regular lon).

    Raises ValueError on any inconsistency. `require_time` should be False for
    pre-computed rate inputs, which have no time dimension.
    """
    required = ("time", "lat", "lon") if require_time else ("lat", "lon")
    for c in required:
        if c not in ds.coords:
            raise ValueError(f"Missing required coordinate: {c}")

    lat_vals = ds["lat"].values
    lon_vals = ds["lon"].values

    if lat_vals.shape[0] != expected_lat_len:
        raise ValueError(
            f"Latitude length {lat_vals.shape[0]} != expected {expected_lat_len}")
    if lon_vals.shape[0] != expected_lon_len:
        raise ValueError(
            f"Longitude length {lon_vals.shape[0]} != expected {expected_lon_len}")

    dlon = np.diff(lon_vals)
    if not np.allclose(dlon, dlon[0], atol=1e-4):
        raise ValueError("Longitude is not regularly spaced (expected regular).")

    dlat = np.diff(lat_vals)
    if np.allclose(dlat, dlat[0], atol=1e-6):
        raise ValueError("Latitude is regularly spaced (expected Gauss-Legendre irregular).")

    print("Dataset coordinates validation passed.")


def convert_units(da: xr.DataArray, from_unit: str, to_unit: str = 'm') -> xr.DataArray:
    """Scale a DataArray between length units (m/dm/cm/mm).

    Divides by (from_scale / to_scale) so values-in-`from_unit` become
    values-in-`to_unit`. Raises ValueError on unknown unit strings.
    """
    if from_unit == to_unit:
        return da

    try:
        from_scale = _UNIT_TO_METERS[from_unit]
        to_scale = _UNIT_TO_METERS[to_unit]
    except KeyError:
        raise ValueError(
            f"Unknown unit; must be one of {sorted(_UNIT_TO_METERS)}. "
            f"Got from={from_unit!r}, to={to_unit!r}.")

    print(f"Converting variable units from {from_unit} to {to_unit}.")
    return da * (from_scale / to_scale)


def standardize_dataset(ds: xr.Dataset,
                        data_var: str,
                        *,
                        input_unit: str = 'm',
                        target_unit: str = 'm',
                        validate: bool = True,
                        aliases: tuple[str, ...] | list[str] | None = None,
                        ) -> xr.Dataset:
    """Full standardization pipeline for a single dataset + variable.

    Applied in order:
      1. `ensure_variable_name` — rename e.g. 'geoid_rate' -> 'delta_g' when
         the canonical name is missing.
      2. `standardize_dims_coords` — lat/lon rename, lon in [0,360], transpose.
      3. `coerce_time_to_years` — normalize the time dtype (skipped for
         pre-computed rate inputs, i.e. no time dim).
      4. Optional `validate_grid` (Gauss-Legendre 257x513).
      5. Optional unit conversion.

    A dataset with no `time` dimension on `data_var` is treated as a
    pre-computed rate: time coercion is skipped, the time-coord check in
    `validate_grid` is relaxed, and the data-variable is tagged with
    `attrs[PRECOMPUTED_RATE_ATTR] = True` so downstream stages (calculate_EWH)
    can skip the rate calculation.
    """
    ds = ensure_variable_name(ds, data_var, aliases=aliases)
    ds = standardize_dims_coords(ds, data_var)

    rate_input = is_precomputed_rate(ds, data_var)
    if rate_input:
        print("WARNING: No time dimension detected. Assuming input grid is "
              "already a rate of change.")
        # netCDF4 attrs can't be Python bool; store as int 1 (round-trips fine).
        ds[data_var].attrs[PRECOMPUTED_RATE_ATTR] = 1
    else:
        ds = coerce_time_to_years(ds)

    if validate:
        validate_grid(ds, require_time=not rate_input)

    if input_unit != target_unit:
        ds[data_var] = convert_units(ds[data_var], input_unit, target_unit)
    return ds


def load_and_standardize(var: str,
                         modelname: str,
                         modeldir: str = '.',
                         *,
                         input_unit: str = 'm',
                         target_unit: str = 'm',
                         validate: bool = True,
                         aliases: tuple[str, ...] | list[str] | None = None,
                         ) -> xr.Dataset:
    """Read `<modeldir>/<modelname>/<var>_<basename>.nc` and standardize it.

    `modelname` may be either a bare model name (joined under `modeldir`) or a
    full path to a per-model directory.
    """
    path = model_input_path(modeldir, modelname, var)
    print(f"Loading {var} for {modelname} from {path}")
    ds = xr.open_dataset(path)
    return standardize_dataset(ds, var,
                               input_unit=input_unit,
                               target_unit=target_unit,
                               validate=validate,
                               aliases=aliases)


def _write_standardized(ds: xr.Dataset,
                        var: str,
                        modelname: str,
                        modeldir: str,
                        outdir_base: str | None,
                        suffix: str,
                        overwrite: bool) -> str:
    """Write standardized dataset next to the raw input (per-model dir).

    If `outdir_base` is None, writes to `<model_dir>/<var>_<basename><suffix>.nc`
    where `<model_dir>` comes from `plotting_utils.model_dir(modeldir, modelname)`.
    If `outdir_base` is given, treats it as a base and writes to
    `<outdir_base>/<basename>/<var>_<basename><suffix>.nc` so the one-file-per-
    model convention is preserved even when the outputs are redirected.
    """
    base = model_basename(modelname)
    if outdir_base is None:
        target_dir = model_dir(modeldir, modelname)
    else:
        target_dir = os.path.join(outdir_base, base)

    os.makedirs(target_dir, exist_ok=True)
    outpath = os.path.join(target_dir, f"{var}_{base}{suffix}.nc")

    if os.path.exists(outpath) and not overwrite:
        print(f"Skipping (exists): {outpath}   (use --overwrite to force)")
        return outpath
    ds.to_netcdf(outpath)
    print(f"Wrote {outpath}")
    return outpath


def run_standardize_stage(variables: Iterable[dict],
                          models: Iterable[str],
                          modeldir: str = '.',
                          outdir: str | None = None,
                          suffix: str = '_std',
                          overwrite: bool = False) -> list[str]:
    """Programmatic entry point used by run_pipeline.py.

    `variables` is an iterable of dicts with keys:
      - name (str, required)
      - input_unit (str, default 'm')
      - target_unit (str, default 'm')
      - validate (bool, default True)
      - aliases (list[str], optional) — extra data-variable names to rename to
        `name` (in addition to DEFAULT_VAR_ALIASES[name]).

    `outdir`, if provided, is treated as a base directory: each model's file
    lands under `<outdir>/<basename>/`. If omitted (default), each standardized
    file is written into its own per-model directory alongside the raw input.
    """
    written: list[str] = []
    for spec in variables:
        var = spec['name']
        input_unit = spec.get('input_unit', 'm')
        target_unit = spec.get('target_unit', 'm')
        validate = bool(spec.get('validate', True))
        aliases = spec.get('aliases')
        for m in models:
            ds = load_and_standardize(var, m, modeldir,
                                      input_unit=input_unit,
                                      target_unit=target_unit,
                                      validate=validate,
                                      aliases=aliases)
            written.append(_write_standardized(ds, var, m, modeldir, outdir,
                                               suffix, overwrite))
    return written


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Standardize <var>_<model>.nc grids to a common convention.")
    p.add_argument("--var", action='append', required=True,
                   help="Variable name (e.g. delta_g). Repeat for multiple variables.")
    p.add_argument("--models", nargs='+', required=True,
                   help="Model names or a text file with one name per line.")
    p.add_argument("--modeldir", type=str, default='.',
                   help=("Parent directory containing per-model subdirectories. "
                         "Each raw input is expected at "
                         "<modeldir>/<model>/<var>_<model>.nc. Bare model names "
                         "are joined under this; full paths are used as-is."))
    p.add_argument("--outdir", type=str, default=None,
                   help=("Optional base directory for standardized outputs. "
                         "If given, each file lands at "
                         "<outdir>/<model>/<var>_<model><suffix>.nc. Defaults "
                         "to writing into the same per-model dir as the input."))
    p.add_argument("--input-unit", type=str, default='m',
                   help="Units of the input variable (m/dm/cm/mm). Default: m.")
    p.add_argument("--target-unit", type=str, default='m',
                   help="Units to write to standardized file. Default: m.")
    p.add_argument("--suffix", type=str, default='_std',
                   help="Suffix appended to <var>_<model> in the output filename.")
    p.add_argument("--overwrite", action='store_true',
                   help="Overwrite existing standardized files.")
    p.add_argument("--no-validate", action='store_true',
                   help="Skip Gauss-Legendre 257x513 grid check.")
    p.add_argument("--alias", action='append', default=None,
                   help=("Additional data-variable name to accept as an alias "
                         "for the canonical --var (repeatable). Applied to "
                         "every --var passed in this invocation."))
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    models = parse_models_arg(args.models)
    variables = [{
        'name': v,
        'input_unit': args.input_unit,
        'target_unit': args.target_unit,
        'validate': not args.no_validate,
        'aliases': args.alias,
    } for v in args.var]

    run_standardize_stage(
        variables=variables,
        models=models,
        modeldir=args.modeldir,
        outdir=args.outdir,
        suffix=args.suffix,
        overwrite=args.overwrite,
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())