#!/usr/bin/env python3
#
# GIAMIP Compliance Checker — check summary
#
# 1. Naming (_check_naming)
#    - Filename parses as {var}_{exp_id}_{group}_{model}.nc, optionally carrying an
#      ensemble member: {var}_{exp_id}_{mYY}_{group}_{model}.nc. The variable name
#      may contain underscores; exp_id, member, group and model must not.
#    - The variable name is a recognised GIAMIP variable name.
#    - The experiment id is one of the recognised GIAMIP experiments
#      (GIAMIP_EXPERIMENTS). The group and model fields are not validated beyond
#      being present and underscore-free: the checker has no roster of groups.
#
# 2. Numerical (_check_numerical)
#    - No NaN or missing values are present in the data array.
#    - No infinite values are present.
#    - No '_FillValue'/'missing_value' attribute is declared (the data request
#      permits no fill value).
#    - All values lie within the variable's plausible "bounds" (loose sanity-check
#      range from the data request; for mask variables this is ~[0, 1]).
#
# 3. Spatial (_check_spatial)  [lat/lon/t variables only]
#    - Grid has exactly 257 latitude × 513 longitude nodes.
#    - Latitude spans [-90, 90] and increases south-to-north.
#    - Longitude spans [0, 360) and increases west-to-east.
#    - The lat/lon values match the canonical GIAMIP grid: Gauss-Legendre
#      latitudes and equally spaced longitudes (see _canonical_lat/_canonical_lon).
#      This needs no forcing file, so it also applies to the experiments whose
#      forcing has not been delivered yet.
#    - If a forcing file is provided, the lat/lon grid matches the forcing grid
#      within a relative tolerance (GIAMIP_GRID_RELTOL).
#
# 3b. Dimension ordering (_check_dimensions)
#    - Variable dimensions are in the expected order: (time, lat, lon) for fields,
#      (time,) for scalars, (degree, order) for Clm/Slm.
#
# 4. Time (_check_time)  [t and lat/lon/t variables only]
#    - The 'time' coordinate variable is present and monotonically increasing.
#    - The companion 'year' variable is present, monotonically increasing, and
#      matches the time axis (year == time/365).
#    - If a forcing file is provided, 'year' is compared step-by-step against the
#      forcing's own 'year' variable.
#
# 5. Attributes (_check_attributes) — reported in three buckets:
#    a) Metadata attributes (global):
#       - group, model, contributors_name, contributors_email present.
#    b) Coordinate attributes:
#       - Time coordinate has units "days since 2000-01-01 00:00:00" and
#         calendar "365_day"; 'year' has units "year".
#       - lat and lon coordinates have units attributes (lat/lon/t variables).
#       - Coordinate variables have the long_name and standard_name requested by
#         the data request (Time/Latitude/Longitude, time/latitude/longitude).
#       - Coordinate variables have the requested data type (float64, except
#         'year' which is int32).
#    c) Variable attributes:
#       - units match the data request.
#       - long_name matches the wording of the data request (case and whitespace
#         insensitive).
#       - standard_name matches the data request (only if one is specified).
#       - dtype matches the precision requested for that variable: float32 for the
#         gridded fields, float64 for the scalar time series and Clm/Slm.
#
# A check that could not be run reports "skipped" and says why; only a check that
# actually ran and found nothing reports "passed".

import argparse
import datetime
import os
import re
import subprocess
import sys

import numpy as np
import xarray as xr
from tqdm import tqdm


GIAMIP_EXPERIMENTS = [f"Exp{i:02d}" for i in range(1, 13)]

# Forcing (ice-history) file associated with each experiment, per experiments.md.
# Several experiments share an ice model but differ in period (122 ka vs 80 ka),
# so the mapping is per-experiment rather than per-ice-model.
# NOTE: only iceHistory-PaleoMIST_1a.nc currently exists; the other names are
# PLACEHOLDERS — update them (or rename the delivered forcing files) to match.
EXPERIMENT_FORCING = {
    "Exp01": "iceHistory-ICE6G_D_122ka.nc",        # ICE-6G_D, 122 ka
    "Exp02": "iceHistory-ICE7G_NA_122ka.nc",       # ICE-7G_NA, 122 ka
    "Exp03": "iceHistory-GLAC3b_profile1_122ka.nc",  # GLAC3b Profile 1, 122 ka
    "Exp04": "iceHistory-GLAC3b_profile2_122ka.nc",  # GLAC3b Profile 2, 122 ka
    "Exp05": "iceHistory-PaleoMIST_1a.nc",         # PaleoMIST (Version a1), 80 ka — available
    "Exp06": "iceHistory-ICE6G_D_80ka.nc",         # ICE-6G_D, 80 ka
    "Exp07": "iceHistory-ICE7G_NA_80ka.nc",        # ICE-7G_NA, 80 ka
    "Exp08": "iceHistory-PaleoMIST_1a.nc",         # PaleoMIST, 80 ka — same as Exp05
    "Exp09": "iceHistory-GLAC3b_profile1_122ka.nc",  # GLAC3b Profile 1, 122 ka — same as Exp03
    "Exp10": "iceHistory-GLAC3b_profile2_122ka.nc",  # GLAC3b Profile 2, 122 ka — same as Exp04
    "Exp11": "iceHistory-GLAC3b_profile3_122ka.nc",  # GLAC3b Profile 3, 122 ka
    "Exp12": "iceHistory-ICE6G_D_122ka.nc",        # ICE-6G_D, 122 ka — same as Exp01
}

# GIAMIP file naming convention:
# {var}_{experiment_id}_{group_name}_{model_name}.nc
GIAMIP_FILENAME_PARTS = 4
GIAMIP_FILENAME_VAR_IDX = 0
GIAMIP_FILENAME_EXP_IDX = 1
GIAMIP_FILENAME_GROUP_IDX = 2
GIAMIP_FILENAME_MODEL_IDX = 3

GIAMIP_GRID_NLAT = 257
GIAMIP_GRID_NLON = 513
GIAMIP_LAT_EXTENT = (-90.0, 90.0)
GIAMIP_LON_EXTENT = (0.0, 360.0)
GIAMIP_COORD_TOL = 0.1  # degrees tolerance for grid extent checks
GIAMIP_GRID_RELTOL = 1e-4  # relative tolerance for matching the forcing lat/lon grid
# Absolute tolerance (degrees) for matching the canonical GIAMIP grid. Independent
# implementations of Gauss-Legendre quadrature agree to ~1e-12 degrees, while any
# grid that is not the prescribed one (equally spaced latitudes, a regional
# subset, a shifted or half-cell-offset grid) is off by ~0.1 degrees or more, so
# 1e-3 degrees (~100 m) separates the two cases with room to spare.
GIAMIP_GRID_ABSTOL = 1e-3


def _canonical_lat() -> np.ndarray:
    """Latitudes of the prescribed GIAMIP grid, south to north.

    The data request prescribes a global Gaussian grid, whose latitudes are the
    Gauss-Legendre quadrature nodes in sin(latitude). The delivered forcing file
    reproduces these to 1e-14 degrees, so they can be regenerated here rather than
    read from a file — which is what lets the grid be validated for the
    experiments whose forcing file has not been delivered yet.
    """
    nodes, _ = np.polynomial.legendre.leggauss(GIAMIP_GRID_NLAT)
    return np.degrees(np.arcsin(nodes))


def _canonical_lon() -> np.ndarray:
    """Longitudes of the prescribed GIAMIP grid, west to east.

    Equally spaced over [0, 360): 513 nodes of 360/513 degrees, so the last node
    is 359.298... and 360 is *not* repeated as 0.
    """
    return np.arange(GIAMIP_GRID_NLON) * (360.0 / GIAMIP_GRID_NLON)


# Time encoding. 'time' is days since 2000-01-01 on a 365-day calendar, and the
# companion 'year' variable holds the year of each step, offset from 2000 CE:
# time 0 -> year 0, 365 -> 1, -365 -> -1.
GIAMIP_YEAR_VAR = "year"
GIAMIP_TIME_CALENDAR = "365_day"
GIAMIP_DAYS_PER_YEAR = 365
# Accepted spellings of the requested time units, "days since 2000-01-01
# 00:00:00". CF treats the shorter forms as identical, so rejecting them would
# fail conformant files.
GIAMIP_TIME_UNITS = "days since 2000-01-01 00:00:00"
_TIME_UNITS_RE = re.compile(
    r"^days\s+since\s+2000-0?1-0?1([\sT]+00:00(:00)?)?\s*$", re.IGNORECASE
)

GIAMIP_VARIABLES = [
    # --- Required 3D (time, lat, lon) ---
    {
        "variable": "delta_bed", "dim": "lat_lon_t", "units": "m",
        "mandatory": True, "output_interval": "forcing",
        "precision": np.float32,
        "standard_name": None, "bounds": (-3000.0, 3000.0),
        "long_name": "Change in the bedrock elevation relative to the initial simulation time step",
    },
    {
        "variable": "delta_g", "dim": "lat_lon_t", "units": "m",
        "mandatory": True, "output_interval": "forcing",
        "precision": np.float32,
        "standard_name": None, "bounds": (-500.0, 500.0),
        "long_name": "Change in the geoid height relative to the initial simulation time step",
    },
    {
        "variable": "ocean_area_fraction", "dim": "lat_lon_t", "units": "1",
        "mandatory": True, "output_interval": "forcing",
        "precision": np.float32,
        "standard_name": "sea_area_fraction", "bounds": (-1e-6, 1.0 + 1e-6),
        "long_name": "Fraction of horizontal grid-cell area covered by ocean",
    },
    {
        "variable": "land_ice_area_fraction", "dim": "lat_lon_t", "units": "1",
        "mandatory": True, "output_interval": "forcing",
        "precision": np.float32,
        "standard_name": "land_ice_area_fraction", "bounds": (-1e-6, 1.0 + 1e-6),
        "long_name": "Fraction of horizontal grid-cell area covered by grounded and floating land ice",
    },
    # --- Required scalars (time only) ---
    {
        "variable": "mean_delta_g", "dim": "t", "units": "m",
        "mandatory": True, "output_interval": "forcing",
        "precision": np.float64,
        "standard_name": None, "bounds": (-300.0, 300.0),
        "long_name": "Spatial mean of geoid height change (delta_g) over the ocean area",
    },
    {
        "variable": "grd_ice_mass", "dim": "t", "units": "kg",
        "mandatory": True, "output_interval": "forcing",
        "precision": np.float64,
        "standard_name": None, "bounds": (1e18, 1e20),
        "long_name": "Spatial integration of grounded ice volume times ice density",
    },
    {
        "variable": "total_ice_mass", "dim": "t", "units": "kg",
        "mandatory": True, "output_interval": "forcing",
        "precision": np.float64,
        "standard_name": None, "bounds": (1e18, 1e20),
        "long_name": "Spatial integration, total (grounded and floating) ice volume times ice density",
    },
    {
        "variable": "ocean_area_grdice", "dim": "t", "units": "m2",
        "mandatory": True, "output_interval": "forcing",
        "precision": np.float64,
        "standard_name": None, "bounds": (1e14, 5e14),
        "long_name": "Total ocean area including marine regions covered by grounded ice",
    },
    {
        "variable": "ocean_area", "dim": "t", "units": "m2",
        "mandatory": True, "output_interval": "forcing",
        "precision": np.float64,
        "standard_name": None, "bounds": (1e14, 5e14),
        "long_name": "Total ocean area excluding marine regions covered by grounded ice",
    },
    {
        "variable": "maf", "dim": "t", "units": "kg",
        "mandatory": True, "output_interval": "forcing",
        "precision": np.float64,
        "standard_name": "land_ice_mass_not_displacing_sea_water", "bounds": (1e19, 1e20),
        "long_name": "Land ice mass above flotation that would contribute to global mean sea-level change if converted to water and added to the ocean",
    },
    # --- Optional 3D ---
    {
        "variable": "delta_bed_east", "dim": "lat_lon_t", "units": "m",
        "mandatory": False, "output_interval": "forcing",
        "precision": np.float32,
        "standard_name": None, "bounds": (-200.0, 200.0),
        "long_name": "Eastward horizontal solid Earth displacement relative to the initial simulation timestep",
    },
    {
        "variable": "delta_bed_north", "dim": "lat_lon_t", "units": "m",
        "mandatory": False, "output_interval": "forcing",
        "precision": np.float32,
        "standard_name": None, "bounds": (-200.0, 200.0),
        "long_name": "Northward horizontal solid Earth displacement relative to the initial simulation timestep",
    },
    # --- Optional spherical harmonics (degree, order) ---
    {
        "variable": "Clm", "dim": "degree_order", "units": "1",
        "mandatory": False, "output_interval": "once",
        "precision": np.float64,
        "standard_name": None, "bounds": None,
        "long_name": "Cosine spherical harmonic coefficients (C_lm) of geoid height change (delta_g) between the first and final simulation timesteps",
    },
    {
        "variable": "Slm", "dim": "degree_order", "units": "1",
        "mandatory": False, "output_interval": "once",
        "precision": np.float64,
        "standard_name": None, "bounds": None,
        "long_name": "Sine spherical harmonic coefficients (S_lm) of geoid height change (delta_g) between the first and final simulation timesteps",
    },
]

GIAMIP_VAR_NAMES = [v["variable"] for v in GIAMIP_VARIABLES]
GIAMIP_MANDATORY_VARS = [v["variable"] for v in GIAMIP_VARIABLES if v["mandatory"]]
# Error count charged to a directory that could not be checked at all (absent,
# empty, or unprocessable): every mandatory variable is missing from it.
_ALL_MANDATORY_MISSING = len(GIAMIP_MANDATORY_VARS)
GIAMIP_VAR_META = {v["variable"]: v for v in GIAMIP_VARIABLES}

# Sort by descending part-count so longer names are tried first (e.g.
# "land_ice_area_fraction" is matched before "land").
_GIAMIP_VAR_NAMES_BY_PARTS = sorted(
    GIAMIP_VAR_NAMES, key=lambda n: n.count("_"), reverse=True
)

# Optional ensemble-member identifier appended to the experiment id, e.g. "m01".
_MEMBER_RE = re.compile(r"m\d+")


def _split_experiment_id(experiment_name: str):
    """Split an experiment id into (base, member).

    "Exp10_m01" -> ("Exp10", "m01");  "Exp05" -> ("Exp05", None).
    """
    base, _, tail = experiment_name.rpartition("_")
    if base and _MEMBER_RE.fullmatch(tail):
        return base, tail
    return experiment_name, None


def _experiment_base(experiment_name: str) -> str:
    """Return the base experiment id, stripping any "_mYY" ensemble member."""
    return _split_experiment_id(experiment_name)[0]


def _parse_giamip_filename(file_name: str):
    """Parse a GIAMIP filename into (var_name, exp_id, group, model) or None.

    Convention: {var}_{exp_id}_{group}_{model}.nc, where exp_id optionally carries
    an ensemble-member identifier: {var}_{exp_id}_{mYY}_{group}_{model}.nc. When a
    member is present it is folded into the returned exp_id (e.g. "Exp10_m01"), so
    each member is treated as a distinct experiment downstream.

    Variable names may contain underscores; exp_id, member, group, and model must
    not. Returns None if no known variable name is found as a prefix, or if the
    trailing fields do not match either the 3-field or member-bearing 4-field form.
    """
    stem = file_name[:-3] if file_name.endswith(".nc") else file_name
    parts = stem.split("_")
    for var in _GIAMIP_VAR_NAMES_BY_PARTS:
        n = var.count("_") + 1  # number of underscore-separated tokens in the var name
        if parts[:n] == var.split("_"):
            trailing = parts[n:]
            if len(trailing) == 3:
                exp, group, model = trailing
                return var, exp, group, model
            if len(trailing) == 4 and _MEMBER_RE.fullmatch(trailing[1]):
                exp, member, group, model = trailing
                return var, f"{exp}_{member}", group, model
    return None


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

EXPECTED_CONDA_ENV = "isschecker"


def warn_if_wrong_environment() -> None:
    """Warn when the script is not running in the environment it was pinned for.

    This is the only environment guard: the checker runs in whichever interpreter
    it is invoked with. Every entry point calls this, so the warning is seen
    however the checker is started.
    """
    active = os.environ.get("CONDA_DEFAULT_ENV", "")
    # Also accept the case where the running interpreter already lives in the
    # expected env (e.g. invoked by absolute path), even if no env is "activated".
    running_in_env = os.path.join("envs", EXPECTED_CONDA_ENV) in sys.executable
    if active != EXPECTED_CONDA_ENV and not running_in_env:
        print(
            f"WARNING: expected conda environment '{EXPECTED_CONDA_ENV}' but"
            f" '{active or '(none)'}' is active. Run 'conda activate {EXPECTED_CONDA_ENV}'"
            " before using this script to ensure the correct package versions are loaded."
        )


def main() -> None:
    warn_if_wrong_environment()
    args = _parse_args()
    run_checker(
        source_path=args.source_dir,
        forcing_path=args.forcing_filepath,
    )


def run_checker(
    source_path: str,
    forcing_path: str | None,
    commit_num: str | None = None,
) -> dict:
    commit_num = _get_commit_number() if commit_num is None else commit_num
    experiments = GIAMIP_EXPERIMENTS
    forcing = _load_forcing(forcing_path) if forcing_path else None

    summary = _run_compliance_checker(
        source_path=source_path,
        commit_num=commit_num,
        experiments=experiments,
        forcing=forcing,
    )

    log_path = os.path.join(source_path, "compliance_checker_log.txt")
    log_text = ""
    if os.path.exists(log_path):
        with open(log_path, "r") as f:
            log_text = f.read()
    summary["log_path"] = log_path
    summary["log_text"] = log_text
    return summary


def _forcing_path_for_experiment(experiment_name: str, forcing_dir: str):
    """Return the forcing file path mapped to an experiment, or None if the
    experiment has no mapping. The returned path is not guaranteed to exist."""
    fname = EXPERIMENT_FORCING.get(_experiment_base(experiment_name))
    if fname is None:
        return None
    return os.path.join(forcing_dir, fname)


def run_model_checker(
    source_dir: str,
    forcing_dir: str,
    commit_num: str | None = None,
) -> dict:
    """Check every experiment subdirectory of a single model directory.

    Each immediate subdirectory whose name is a recognised experiment (e.g.
    'Exp05' or 'Exp08_m01') is checked with the forcing file mapped to that
    experiment in EXPERIMENT_FORCING, looked up inside forcing_dir. If the
    forcing file is missing, the experiment is still checked but the
    forcing-based time/grid comparisons are skipped.

    Returns a dict mapping experiment id -> summary.
    """
    commit_num = _get_commit_number() if commit_num is None else commit_num

    if not os.path.isdir(source_dir):
        print(f"ERROR: Model directory not found: '{source_dir}'.")
        return {}

    results = {}
    for name in sorted(os.listdir(source_dir)):
        exp_dir = os.path.join(source_dir, name)
        if not os.path.isdir(exp_dir):
            continue
        if _experiment_base(name) not in GIAMIP_EXPERIMENTS:
            print(f"Skipping '{name}': not a recognised GIAMIP experiment directory.")
            continue

        forcing_path = _forcing_path_for_experiment(name, forcing_dir)
        if forcing_path is None:
            print(f"WARNING: no forcing mapping for experiment '{name}'.")
        elif not os.path.exists(forcing_path):
            print(
                f"WARNING: forcing file '{forcing_path}' for experiment '{name}'"
                " not found; running without forcing-based time/grid checks."
            )
            forcing_path = None

        results[name] = run_checker(
            source_path=exp_dir, forcing_path=forcing_path,
            commit_num=commit_num,
        )

    _print_rollup(source_dir, results)
    return results


def run_full_checker(
    output_dir: str,
    forcing_dir: str,
    commit_num: str | None = None,
) -> dict:
    """Check every model directory under output_dir, each looping over its
    experiments (see run_model_checker). Returns a dict mapping model name ->
    {experiment id: summary}.
    """
    commit_num = _get_commit_number() if commit_num is None else commit_num

    if not os.path.isdir(output_dir):
        print(f"ERROR: Output directory not found: '{output_dir}'.")
        return {}

    results = {}
    for model in sorted(os.listdir(output_dir)):
        model_dir = os.path.join(output_dir, model)
        if not os.path.isdir(model_dir):
            continue
        print(f"\n######## Model: {model} ########")
        results[model] = run_model_checker(
            model_dir, forcing_dir, commit_num=commit_num,
        )
    return results


def _print_rollup(source_dir: str, results: dict) -> None:
    print("\n=========================================================================")
    print(f"Roll-up for {source_dir}:")
    if not results:
        print("  (no experiment directories checked)")
    for exp in sorted(results):
        print(f"  {exp:16} : {results[exp]['total_errors']} error(s)")
    print("=========================================================================")


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_commit_number() -> str:
    try:
        process = subprocess.Popen(
            ["git", "log", "--pretty=format:%h", "-n", "1"], stdout=subprocess.PIPE
        )
        commit_num, _ = process.communicate()
        return commit_num.decode("UTF-8")
    except Exception:
        return "No commit number identified."


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check GIAMIP simulation NetCDF datasets for compliance."
    )
    parser.add_argument(
        "--source-dir",
        required=True,
        help="Directory containing the GIAMIP NetCDF files to check.",
    )
    parser.add_argument(
        "--forcing-filepath",
        required=True,
        help="Path to the forcing NetCDF file used to validate output time axes.",
    )
    return parser.parse_args()


def _load_forcing(forcing_path: str):
    """Load the forcing file's reference time vector and spatial grid.

    Returns a dict with keys:
      - "years": integer array of calendar years (from the 'year' variable), or None.
      - "lat":   forcing latitude array (float), or None.
      - "lon":   forcing longitude array (float), or None.

    The forcing time axis is encoded with a problematic calendar (e.g. "days since
    2000" / "365_day"), so we use the explicit 'year' variable as the reference time
    vector instead of decoding 'time'. Returns None if the file cannot be opened.
    """
    try:
        ds = xr.open_dataset(forcing_path, decode_times=False)
    except Exception as e:
        print(f"WARNING: could not load forcing file '{forcing_path}': {e}")
        return None

    if "year" in ds.variables:
        years = np.asarray(ds["year"].values).astype(int)
    else:
        print(f"WARNING: forcing file '{forcing_path}' has no 'year' variable;"
              " time axis comparison will be skipped.")
        years = None

    if "lat" in ds.variables and "lon" in ds.variables:
        lat = np.asarray(ds["lat"].values, dtype=float)
        lon = np.asarray(ds["lon"].values, dtype=float)
    else:
        print(f"WARNING: forcing file '{forcing_path}' has no 'lat'/'lon' grid;"
              " spatial grid comparison will be skipped.")
        lat = lon = None

    return {"years": years, "lat": lat, "lon": lon}


def _run_compliance_checker(
    source_path: str,
    commit_num: str,
    experiments: list,
    forcing=None,
) -> dict:
    if not os.path.isdir(source_path):
        print(f"ERROR: Directory not found: '{source_path}'.")
        return _empty_summary(_ALL_MANDATORY_MISSING)

    try:
        log_path = os.path.join(source_path, "compliance_checker_log.txt")
        with open(log_path, "w") as f:
            print("-> Checking " + source_path)
            print()
            _write_log_header(f, commit_num, source_path, datetime.date.today())

            experiment_groups = _group_files_by_experiment(source_path)
            if not experiment_groups:
                # An empty directory is not a compliant one. Report every mandatory
                # variable as missing so that an experiment which was never
                # delivered — or whose transfer silently failed — can never come
                # back as "successfully verified with no errors".
                msg = (
                    f"No .nc files found in directory '{source_path}'."
                    f" All {_ALL_MANDATORY_MISSING} mandatory variable(s) are"
                    f" therefore missing: {GIAMIP_MANDATORY_VARS}"
                )
                print(f"ERROR: {msg}")
                f.write(f"ERROR: {msg}\n")
                summary = _empty_summary(_ALL_MANDATORY_MISSING)
                _print_total_summary(source_path, summary["total_errors"])
            else:
                summary = _process_experiments(
                    log_file=f,
                    source_path=source_path,
                    experiment_groups=experiment_groups,
                    experiments=experiments,
                    forcing=forcing,
                )

        _insert_synthesis(
            source_path=source_path,
            exp_counter=summary["exp_counter"],
            file_counter=summary["file_counter"],
            total_errors=summary["total_errors"],
            total_file_errors=summary["total_file_errors"],
            total_naming_errors=summary["total_naming_errors"],
            total_num_errors=summary["total_num_errors"],
            total_spatial_errors=summary["total_spatial_errors"],
            total_dim_errors=summary["total_dim_errors"],
            total_time_errors=summary["total_time_errors"],
            total_meta_attr_errors=summary["total_meta_attr_errors"],
            total_coord_attr_errors=summary["total_coord_attr_errors"],
            total_var_attr_errors=summary["total_var_attr_errors"],
            report_naming_issues=summary["report_naming_issues"],
        )
        return summary

    except TypeError as err:
        print("Something went wrong. Error:", err)
        return _empty_summary(_ALL_MANDATORY_MISSING)


def _empty_summary(file_errors: int = 0) -> dict:
    """Summary for a directory that could not be checked file by file.

    *file_errors* is attributed to the "missing mandatory variable" bucket: a
    directory that is absent, empty, or that the checker could not process is
    treated as one in which every mandatory variable is missing. Reporting zero
    errors there would be a false pass — the worst failure mode for this tool —
    so the count is never silently zero.
    """
    return {
        "exp_counter": 0,
        "file_counter": 0,
        "total_errors": file_errors,
        "total_naming_errors": 0,
        "total_num_errors": 0,
        "total_spatial_errors": 0,
        "total_dim_errors": 0,
        "total_time_errors": 0,
        "total_meta_attr_errors": 0,
        "total_coord_attr_errors": 0,
        "total_var_attr_errors": 0,
        "total_file_errors": file_errors,
        "report_naming_issues": [],
    }


def _group_files_by_experiment(source_path: str) -> dict:
    groups = {}
    for fname in sorted(os.listdir(source_path)):
        if not fname.endswith(".nc"):
            continue
        parsed = _parse_giamip_filename(fname)
        exp_name = parsed[1] if parsed is not None else "_unknown"
        groups.setdefault(exp_name, []).append(fname)
    return groups


def _process_experiments(
    log_file,
    source_path: str,
    experiment_groups: dict,
    experiments: list,
    forcing=None,
) -> dict:
    total_naming_errors = 0
    total_num_errors = 0
    total_spatial_errors = 0
    total_dim_errors = 0
    total_time_errors = 0
    total_meta_attr_errors = 0
    total_coord_attr_errors = 0
    total_var_attr_errors = 0
    total_file_errors = 0
    report_naming_issues = []
    file_counter = 0
    exp_counter = 0

    for experiment_name, exp_files in experiment_groups.items():
        exp_counter += 1
        exp_summary = _process_single_experiment(
            log_file=log_file,
            source_path=source_path,
            experiment_name=experiment_name,
            exp_files=exp_files,
            experiments=experiments,
            report_naming_issues=report_naming_issues,
            forcing=forcing,
        )
        file_counter += exp_summary["file_counter"]
        total_naming_errors += exp_summary["exp_naming_errors"]
        total_num_errors += exp_summary["exp_num_errors"]
        total_spatial_errors += exp_summary["exp_spatial_errors"]
        total_dim_errors += exp_summary["exp_dim_errors"]
        total_time_errors += exp_summary["exp_time_errors"]
        total_meta_attr_errors += exp_summary["exp_meta_attr_errors"]
        total_coord_attr_errors += exp_summary["exp_coord_attr_errors"]
        total_var_attr_errors += exp_summary["exp_var_attr_errors"]
        total_file_errors += exp_summary["exp_file_errors"]

        _print_experiment_summary(exp_summary["experiment_name"], exp_summary["exp_errors"])

    total_errors = (
        total_naming_errors + total_num_errors + total_spatial_errors + total_dim_errors
        + total_time_errors + total_meta_attr_errors + total_coord_attr_errors
        + total_var_attr_errors + total_file_errors
    )
    _print_total_summary(source_path, total_errors)

    return {
        "exp_counter": exp_counter,
        "file_counter": file_counter,
        "total_errors": total_errors,
        "total_naming_errors": total_naming_errors,
        "total_num_errors": total_num_errors,
        "total_spatial_errors": total_spatial_errors,
        "total_dim_errors": total_dim_errors,
        "total_time_errors": total_time_errors,
        "total_meta_attr_errors": total_meta_attr_errors,
        "total_coord_attr_errors": total_coord_attr_errors,
        "total_var_attr_errors": total_var_attr_errors,
        "total_file_errors": total_file_errors,
        "report_naming_issues": report_naming_issues,
    }


def _process_single_experiment(
    log_file,
    source_path: str,
    experiment_name: str,
    exp_files: list,
    experiments: list,
    report_naming_issues: list,
    forcing=None,
) -> dict:
    exp_naming_errors = 0
    exp_num_errors = 0
    exp_spatial_errors = 0
    exp_dim_errors = 0
    exp_time_errors = 0
    exp_meta_attr_errors = 0
    exp_coord_attr_errors = 0
    exp_var_attr_errors = 0
    exp_file_errors = 0

    log_file.write("\n ")
    log_file.write("**********************************************************\n")
    log_file.write(f" ** Experiment: {experiment_name} \n ")
    log_file.write("**********************************************************\n")
    log_file.write("\n ")

    if _experiment_base(experiment_name) not in experiments:
        log_file.write(
            f"ERROR: The compliance check is ignored for experiment '{experiment_name}'"
            f" as it is not in {experiments}.\n"
        )
        exp_naming_errors += 1
        report_naming_issues.append(
            f"Compliance check ignored: experiment '{experiment_name}' not in the experiments list."
        )
        return {
            "file_counter": 0,
            "experiment_name": experiment_name,
            "exp_errors": exp_naming_errors,
            "exp_naming_errors": exp_naming_errors,
            "exp_num_errors": 0,
            "exp_spatial_errors": 0,
            "exp_dim_errors": 0,
            "exp_time_errors": 0,
            "exp_meta_attr_errors": 0,
            "exp_coord_attr_errors": 0,
            "exp_var_attr_errors": 0,
            "exp_file_errors": 0,
        }

    # Check mandatory variables
    present_vars = set()
    for f in exp_files:
        parsed = _parse_giamip_filename(f)
        if parsed is not None:
            present_vars.add(parsed[0])
    missing_mandatory = [v for v in GIAMIP_MANDATORY_VARS if v not in present_vars]
    if not missing_mandatory:
        log_file.write(
            f"Missing mandatory variable test: {experiment_name}: all mandatory variables present.\n"
        )
    else:
        log_file.write(
            f"ERROR: In experiment {experiment_name}, these mandatory variable(s) are missing:"
            f" {missing_mandatory}\n"
        )
        exp_file_errors += len(missing_mandatory)

    file_counter = 0
    for fname in tqdm(exp_files):
        file_counter += 1
        file_summary = _process_single_file(
            log_file=log_file,
            source_path=source_path,
            file=fname,
            experiment_name=experiment_name,
            experiments=experiments,
            report_naming_issues=report_naming_issues,
            forcing=forcing,
        )
        exp_naming_errors += file_summary["var_naming_errors"]
        exp_num_errors += file_summary["var_num_errors"]
        exp_spatial_errors += file_summary["var_spatial_errors"]
        exp_dim_errors += file_summary["var_dim_errors"]
        exp_time_errors += file_summary["var_time_errors"]
        exp_meta_attr_errors += file_summary["var_meta_attr_errors"]
        exp_coord_attr_errors += file_summary["var_coord_attr_errors"]
        exp_var_attr_errors += file_summary["var_var_attr_errors"]

    exp_errors = (
        exp_naming_errors + exp_num_errors + exp_spatial_errors + exp_dim_errors
        + exp_time_errors + exp_meta_attr_errors + exp_coord_attr_errors
        + exp_var_attr_errors + exp_file_errors
    )
    return {
        "file_counter": file_counter,
        "experiment_name": experiment_name,
        "exp_errors": exp_errors,
        "exp_naming_errors": exp_naming_errors,
        "exp_num_errors": exp_num_errors,
        "exp_spatial_errors": exp_spatial_errors,
        "exp_dim_errors": exp_dim_errors,
        "exp_time_errors": exp_time_errors,
        "exp_meta_attr_errors": exp_meta_attr_errors,
        "exp_coord_attr_errors": exp_coord_attr_errors,
        "exp_var_attr_errors": exp_var_attr_errors,
        "exp_file_errors": exp_file_errors,
    }


def _process_single_file(
    log_file,
    source_path: str,
    file: str,
    experiment_name: str,
    experiments: list,
    report_naming_issues: list,
    forcing=None,
) -> dict:
    zero = {"var_naming_errors": 0, "var_num_errors": 0,
            "var_spatial_errors": 0, "var_dim_errors": 0, "var_time_errors": 0,
            "var_meta_attr_errors": 0, "var_coord_attr_errors": 0,
            "var_var_attr_errors": 0}

    parsed = _parse_giamip_filename(file)
    var_name = parsed[0] if parsed is not None else file[:-3].split("_")[0]

    # Open dataset (decode_times=False so we can inspect encoding)
    try:
        ds = xr.open_dataset(os.path.join(source_path, file), decode_times=False)
    except Exception as e:
        log_file.write(f" - ERROR: Cannot open {file}: {e}\n")
        return {**zero, "var_naming_errors": 1}

    (naming_errors, num_errors, spatial_errors, dim_errors, time_errors,
     meta_attr_errors, coord_attr_errors, var_attr_errors) = _run_variable_checks(
        log_file=log_file,
        ds=ds,
        file_name=file,
        var_name=var_name,
        experiment_name=experiment_name,
        experiments=experiments,
        report_naming_issues=report_naming_issues,
        forcing=forcing,
    )

    total = (naming_errors + num_errors + spatial_errors + dim_errors + time_errors
             + meta_attr_errors + coord_attr_errors + var_attr_errors)
    log_file.write("\n")
    log_file.write("----------------------------------------------------------\n")
    log_file.write(f"{experiment_name} - {var_name} - File: {file}\n")
    if total > 0:
        log_file.write(f"{total} error(s). Please review before sharing.\n")
    else:
        log_file.write("No errors. Good job !\n")
    log_file.write("No warnings.\n")
    log_file.write("----------------------------------------------------------\n")

    return {
        "var_naming_errors": naming_errors,
        "var_num_errors": num_errors,
        "var_spatial_errors": spatial_errors,
        "var_dim_errors": dim_errors,
        "var_time_errors": time_errors,
        "var_meta_attr_errors": meta_attr_errors,
        "var_coord_attr_errors": coord_attr_errors,
        "var_var_attr_errors": var_attr_errors,
    }


# ---------------------------------------------------------------------------
# Check functions
# ---------------------------------------------------------------------------

def _check_naming(
    log_file,
    file_name: str,
    var_name: str,
    experiment_name: str,
    experiments: list,
    report_naming_issues: list,
) -> int:
    errors = 0
    log_file.write("NAMING Tests \n")

    parsed = _parse_giamip_filename(file_name)
    if parsed is None:
        log_file.write(
            f" - ERROR: filename '{file_name}' does not follow the convention"
            f" {{var}}_{{exp_id}}_{{group}}_{{model}}.nc — no known GIAMIP variable"
            f" found as a prefix, or wrong number of trailing fields.\n"
        )
        report_naming_issues.append(f"Compliance check ignored: '{file_name}' wrong field count.")
        return errors + 1

    if var_name not in GIAMIP_VAR_NAMES:
        log_file.write(
            f" - ERROR: variable '{var_name}' is not a recognised GIAMIP variable name.\n"
        )
        errors += 1

    if _experiment_base(experiment_name) not in experiments:
        log_file.write(
            f" - ERROR: experiment '{experiment_name}' is not a known GIAMIP experiment.\n"
        )
        errors += 1

    if errors == 0:
        log_file.write(f" - Filename convention: passed\n")
    return errors


# Attributes that declare a missing/invalid value. The data request permits none
# of them: "No invalid or missing values are permitted for this variable."
_FILL_VALUE_ATTRS = ("_FillValue", "missing_value")


def _check_no_fill_value(log_file, ds: xr.Dataset, var_name: str) -> int:
    """Report any missing-value attribute declared on *var_name*.

    xarray moves '_FillValue' from the file into '.encoding' when it decodes the
    array, so both encoding and attrs are inspected.
    """
    errors = 0
    declared = {**ds[var_name].encoding, **ds[var_name].attrs}
    for attr in _FILL_VALUE_ATTRS:
        value = declared.get(attr)
        if value is not None:
            log_file.write(
                f" - ERROR: variable '{var_name}' declares '{attr}' = {value};"
                f" the data request permits no fill value.\n"
            )
            errors += 1
    if errors == 0:
        log_file.write(" - No fill value declared: passed\n")
    return errors


def _check_numerical(
    log_file,
    ds: xr.Dataset,
    var_name: str,
) -> int:
    errors = 0
    log_file.write("NUMERICAL Tests \n")

    if var_name not in ds:
        log_file.write(f" - ERROR: variable '{var_name}' not found in dataset.\n")
        return errors + 1

    # No NaN / missing values
    data = ds[var_name].values
    if np.isnan(data).any():
        log_file.write(f" - ERROR: variable '{var_name}' contains NaN / missing values.\n")
        errors += 1
    else:
        log_file.write(f" - No missing values: passed\n")

    # No infinities (the data request permits no invalid values either).
    if np.isinf(data).any():
        log_file.write(f" - ERROR: variable '{var_name}' contains infinite values.\n")
        errors += 1

    # No fill value may be declared: the data request specifies _FillValue "none"
    # for every variable. Values matching a declared _FillValue have already been
    # decoded to NaN by xarray and are caught above; this catches the declaration
    # itself, which signals that missing data is expected.
    errors += _check_no_fill_value(log_file, ds, var_name)

    # Values must lie within the variable's plausible range (sanity bounds).
    bounds = GIAMIP_VAR_META.get(var_name, {}).get("bounds")
    if bounds is not None:
        lo, hi = bounds
        vmin = float(np.nanmin(data))
        vmax = float(np.nanmax(data))
        if vmin < lo or vmax > hi:
            log_file.write(
                f" - ERROR: variable '{var_name}' has values outside its plausible"
                f" range [{lo}, {hi}] (min={vmin}, max={vmax}).\n"
            )
            errors += 1
        else:
            log_file.write(f" - Value range within [{lo}, {hi}]: passed\n")

    return errors


def _check_spatial(log_file, ds: xr.Dataset, forcing_lat=None, forcing_lon=None) -> int:
    errors = 0
    log_file.write("SPATIAL Tests \n")

    for coord in ("lat", "lon"):
        if coord not in ds.coords:
            log_file.write(f" - ERROR: coordinate '{coord}' not found.\n")
            errors += 1

    if errors:
        return errors

    nlat = int(ds["lat"].size)
    nlon = int(ds["lon"].size)

    if nlat == GIAMIP_GRID_NLAT:
        log_file.write(f" - Latitude size ({nlat} nodes): passed\n")
    else:
        log_file.write(
            f" - ERROR: latitude has {nlat} nodes, expected {GIAMIP_GRID_NLAT}.\n"
        )
        errors += 1

    if nlon == GIAMIP_GRID_NLON:
        log_file.write(f" - Longitude size ({nlon} nodes): passed\n")
    else:
        log_file.write(
            f" - ERROR: longitude has {nlon} nodes, expected {GIAMIP_GRID_NLON}.\n"
        )
        errors += 1

    lat_vals = ds["lat"].values.astype(float)
    lon_vals = ds["lon"].values.astype(float)

    lat_min, lat_max = float(lat_vals.min()), float(lat_vals.max())
    lon_min, lon_max = float(lon_vals.min()), float(lon_vals.max())

    if lat_min > GIAMIP_LAT_EXTENT[0]:
        log_file.write(f" - Latitude south edge ({lat_min}° > -90°): passed\n")
    else:
        log_file.write(
            f" - ERROR: latitude south edge {lat_min}°, expected > {GIAMIP_LAT_EXTENT[0]}°.\n"
        )
        errors += 1

    if lat_max < GIAMIP_LAT_EXTENT[1]:
        log_file.write(f" - Latitude north edge ({lat_max}° < 90°): passed\n")
    else:
        log_file.write(
            f" - ERROR: latitude north edge {lat_max}°, expected < {GIAMIP_LAT_EXTENT[1]}°.\n"
        )
        errors += 1

    if abs(lon_min - GIAMIP_LON_EXTENT[0]) <= GIAMIP_COORD_TOL:
        log_file.write(f" - Longitude west edge ({lon_min}°): passed\n")
    else:
        log_file.write(
            f" - ERROR: longitude west edge {lon_min}°, expected {GIAMIP_LON_EXTENT[0]}°.\n"
        )
        errors += 1

    if lon_max < GIAMIP_LON_EXTENT[1]:
        log_file.write(f" - Longitude east edge ({lon_max}° < 360°): passed\n")
    else:
        log_file.write(
            f" - ERROR: longitude east edge {lon_max}° must be < {GIAMIP_LON_EXTENT[1]}°.\n"
        )
        errors += 1

    if _strictly_increasing(lat_vals):
        log_file.write(" - Latitude increases south-to-north: passed\n")
    else:
        log_file.write(" - ERROR: latitude is not monotonically increasing (south-to-north).\n")
        errors += 1

    if _strictly_increasing(lon_vals):
        log_file.write(" - Longitude increases west-to-east: passed\n")
    else:
        log_file.write(" - ERROR: longitude is not monotonically increasing (west-to-east).\n")
        errors += 1

    # Grid values must be the prescribed GIAMIP grid. This is checked against the
    # analytically regenerated grid, so it does not depend on a forcing file being
    # available; the forcing comparison below is the same check against the file
    # actually used to drive the experiment.
    errors += _check_matches_canonical_grid(log_file, "latitude", lat_vals, _canonical_lat())
    errors += _check_matches_canonical_grid(log_file, "longitude", lon_vals, _canonical_lon())

    # Grid must match the forcing/input grid (within a relative tolerance).
    errors += _check_grid_matches_forcing(log_file, "latitude", lat_vals, forcing_lat)
    errors += _check_grid_matches_forcing(log_file, "longitude", lon_vals, forcing_lon)

    return errors


def _check_matches_canonical_grid(log_file, name: str, vals, canonical) -> int:
    """Verify one coordinate axis against the prescribed GIAMIP grid.

    The node-count and extent checks above pass for grids that are the right size
    but the wrong grid — a regional subset, an equally spaced (non-Gaussian)
    latitude axis, or a half-cell-shifted longitude axis all clear them. This
    compares the values themselves. Returns 0 (skipped) if the axis has the wrong
    number of nodes, since that is already reported.
    """
    vals = np.asarray(vals, dtype=float)
    if vals.shape != canonical.shape:
        log_file.write(
            f" - {name.capitalize()} grid vs the GIAMIP grid: skipped"
            f" (wrong number of nodes, reported above).\n"
        )
        return 0

    max_diff = float(np.max(np.abs(vals - canonical)))
    if max_diff <= GIAMIP_GRID_ABSTOL:
        log_file.write(f" - {name.capitalize()} matches the GIAMIP grid: passed\n")
        return 0

    log_file.write(
        f" - ERROR: {name} does not match the prescribed GIAMIP grid"
        f" (max abs difference {max_diff:.6g}°, tolerance {GIAMIP_GRID_ABSTOL}°).\n"
        f"   Expected {canonical[0]:.6f}° ... {canonical[-1]:.6f}°,"
        f" found {vals[0]:.6f}° ... {vals[-1]:.6f}°.\n"
    )
    _write_grid_hint(log_file, name, vals, canonical)
    return 1


def _write_grid_hint(log_file, name: str, vals, canonical) -> None:
    """Name the likely cause of a canonical-grid mismatch, when it is recognisable."""
    spacing = np.diff(vals)
    covered = float(vals[-1] - vals[0])
    expected_cover = float(canonical[-1] - canonical[0])

    if covered < 0.9 * expected_cover:
        log_file.write(
            f"   HINT: the {name} axis covers only {covered:.3f}° of the"
            f" {expected_cover:.3f}° expected — this looks like a regional grid, or a\n"
            "   global field that was cropped. GIAMIP output must be on the full global grid.\n"
        )
    elif name == "latitude" and np.allclose(spacing, spacing[0], atol=1e-6):
        log_file.write(
            "   HINT: the latitude axis is equally spaced. The GIAMIP grid is Gaussian:\n"
            "   its spacing varies from 0.693° near the poles to 0.699° near the equator.\n"
            "   Interpolate onto the forcing file's lat/lon axes rather than building an\n"
            "   equally spaced axis with linspace/arange.\n"
        )


def _check_grid_matches_forcing(log_file, name: str, vals, forcing_vals) -> int:
    """Verify a coordinate axis matches the forcing grid within GIAMIP_GRID_RELTOL.

    A small absolute tolerance (scaled by the relative tolerance) is included so that
    coordinates passing through zero are compared sensibly. Returns 0 if no forcing
    grid is available (comparison skipped).
    """
    if forcing_vals is None:
        return 0
    forcing_vals = np.asarray(forcing_vals, dtype=float)
    if vals.shape != forcing_vals.shape:
        log_file.write(
            f" - ERROR: {name} grid has {vals.shape[0]} node(s),"
            f" forcing grid has {forcing_vals.shape[0]} node(s).\n"
        )
        return 1
    atol = GIAMIP_GRID_RELTOL * float(np.max(np.abs(forcing_vals)))
    if np.allclose(vals, forcing_vals, rtol=GIAMIP_GRID_RELTOL, atol=atol):
        log_file.write(f" - {name.capitalize()} grid matches forcing: passed\n")
        return 0
    max_diff = float(np.max(np.abs(vals - forcing_vals)))
    log_file.write(
        f" - ERROR: {name} grid does not match the forcing grid"
        f" (max abs difference {max_diff}, reltol {GIAMIP_GRID_RELTOL}).\n"
    )
    return 1


# Shown when the dimensions are in exactly the reverse of the expected order,
# which is the signature of writing the file from a column-major language: the
# netCDF Fortran and MATLAB interfaces reverse the dimension list relative to the
# C/CDL view, so an order transcribed literally from the data request comes out
# reversed in the file itself.
_COLUMN_MAJOR_HINT = (
    "   HINT: the dimensions are in exactly the reverse of the expected order.\n"
    "   If you write NetCDF from a column-major language (Fortran, MATLAB), those\n"
    "   interfaces reverse the dimension list relative to the C/Python view that the\n"
    "   data request is written in, so declaring\n"
    "       {expected}\n"
    "   there produces\n"
    "       {reversed}\n"
    "   in the file itself. Declare them in the reverse order in Fortran/MATLAB, or\n"
    "   transpose the array before writing. Either way, confirm with\n"
    "   'ncdump -h <file>': it shows the order as this checker sees it.\n"
)

# Expected dimension ordering for each dim type. The time dimension may be named
# "time" or "t" in the data; both are normalised to "time" before comparison.
_EXPECTED_DIMS = {
    "lat_lon_t": ("time", "lat", "lon"),
    "t": ("time",),
    "degree_order": ("degree", "order"),
}


def _check_dimensions(log_file, ds: xr.Dataset, var_name: str, dim_type: str) -> int:
    errors = 0
    log_file.write("DIMENSION Tests \n")

    expected = _EXPECTED_DIMS.get(dim_type)
    if var_name not in ds:
        log_file.write(
            f" - Dimension ordering: skipped (variable '{var_name}' is not in the file).\n"
        )
        return errors
    if expected is None:
        log_file.write(
            f" - Dimension ordering: skipped (no expected order defined for"
            f" dimension type '{dim_type}').\n"
        )
        return errors

    actual = tuple("time" if d == "t" else d for d in ds[var_name].dims)
    if actual == expected:
        log_file.write(f" - Dimension ordering {expected}: passed\n")
    else:
        log_file.write(
            f" - ERROR: dimension ordering {tuple(ds[var_name].dims)}"
            f" does not match expected {expected}.\n"
        )
        errors += 1
        if actual == tuple(reversed(expected)):
            log_file.write(_COLUMN_MAJOR_HINT.format(expected=expected,
                                                     reversed=tuple(reversed(expected))))
    return errors


# Accepted names for the time axis. "time" is what the data request asks for;
# "t" is accepted so that a file using it is checked rather than skipped, and is
# reported by the attribute checks.
_TIME_COORD_NAMES = ("time", "t")


def _find_time_coord(ds: xr.Dataset):
    """Return the name of the time coordinate *variable*, or None.

    Looked up in ds.variables rather than ds.dims on purpose: when a dimension has
    no coordinate variable, xarray fabricates an index (0, 1, 2, ...) for it, so
    ds["time"] would silently return a made-up axis and every time check would
    "pass" on a file that carries no time axis at all.
    """
    return next((n for n in _TIME_COORD_NAMES if n in ds.variables), None)


def _check_year_matches_time(log_file, time_vals, var_years) -> int:
    """Check that 'year' is the year of each 'time' step.

    'time' is days since 2000-01-01 on a 365-day calendar, so the year of a step
    is time/365 (0 -> 0, 365 -> 1, -365 -> -1). 'year' is an integer, so it is
    accepted as the nearest year: half a year of slack. The usual failure this
    catches is years written into 'time' itself.
    """
    if len(time_vals) != len(var_years):
        log_file.write(
            f" - ERROR: '{GIAMIP_YEAR_VAR}' has {len(var_years)} step(s),"
            f" 'time' has {len(time_vals)} step(s).\n"
        )
        return 1

    implied = time_vals / GIAMIP_DAYS_PER_YEAR
    deviation = np.abs(implied - var_years)
    worst = int(np.argmax(deviation))
    if deviation[worst] > 0.5:
        log_file.write(
            f" - ERROR: '{GIAMIP_YEAR_VAR}' does not match 'time'/{GIAMIP_DAYS_PER_YEAR}"
            f" (worst at index {worst}: time {time_vals[worst]} implies year"
            f" {implied[worst]:.3f}, but {GIAMIP_YEAR_VAR} is {var_years[worst]}).\n"
        )
        return 1

    log_file.write(f" - '{GIAMIP_YEAR_VAR}' matches the 'time' axis: passed\n")
    return 0


def _check_time(
    log_file,
    ds: xr.Dataset,
    output_interval: str,
    forcing_years=None,
) -> int:
    errors = 0
    log_file.write("TIME Tests \n")

    time_coord = _find_time_coord(ds)
    if time_coord is None:
        if any(d in ds.dims for d in _TIME_COORD_NAMES):
            log_file.write(
                " - ERROR: the file has a time dimension but no 'time' coordinate"
                " variable; the time axis itself is missing.\n"
            )
        else:
            log_file.write(" - ERROR: no time dimension and no 'time' coordinate variable.\n")
        return errors + 1

    time_vals = np.asarray(ds[time_coord].values, dtype=float)
    if not _strictly_increasing(time_vals):
        log_file.write(" - ERROR: time is not monotonically increasing.\n")
        errors += 1
    else:
        log_file.write(" - Time is monotonically increasing: passed\n")

    # The reference time vector is the 'year' variable, not the 'time' axis:
    # 'time' is days since 2000-01-01 on a 365-day calendar, so year == time/365.
    if GIAMIP_YEAR_VAR not in ds.variables:
        log_file.write(
            f" - ERROR: variable '{GIAMIP_YEAR_VAR}' not found; it is required"
            f" alongside 'time' and carries the year of each output step.\n"
        )
        return errors + 1

    try:
        var_years = np.asarray(ds[GIAMIP_YEAR_VAR].values).astype(int)
    except Exception as err:
        log_file.write(
            f" - ERROR: '{GIAMIP_YEAR_VAR}' could not be read as years: {err}\n"
        )
        return errors + 1

    if not _strictly_increasing(var_years):
        log_file.write(f" - ERROR: '{GIAMIP_YEAR_VAR}' is not monotonically increasing.\n")
        return errors + 1
    log_file.write(f" - '{GIAMIP_YEAR_VAR}' is monotonically increasing: passed\n")

    errors += _check_year_matches_time(log_file, time_vals, var_years)

    # Compare against the forcing 'year' variable, year-by-year.
    if output_interval == "forcing" and forcing_years is not None:
        if len(var_years) != len(forcing_years):
            log_file.write(
                f" - ERROR: '{GIAMIP_YEAR_VAR}' has {len(var_years)} step(s),"
                f" forcing has {len(forcing_years)} step(s).\n"
            )
            errors += 1
        elif not np.array_equal(var_years, forcing_years):
            diffs = np.where(var_years != forcing_years)[0]
            log_file.write(
                f" - ERROR: time axis does not match forcing at {len(diffs)} step(s)"
                f" (first mismatch: index {diffs[0]},"
                f" variable year {var_years[diffs[0]]},"
                f" forcing year {forcing_years[diffs[0]]}).\n"
            )
            errors += 1
        else:
            log_file.write(" - Time axis matches forcing: passed\n")

    return errors


# CF standard_name requested for each coordinate variable by the data request.
GIAMIP_COORD_STANDARD_NAMES = {
    "time": "time",
    "lat": "latitude",
    "lon": "longitude",
}

# long_name requested for each coordinate variable by the data request.
GIAMIP_COORD_LONG_NAMES = {
    "time": "Time",
    "year": "Calendar year corresponding to each output time"
            " (positive = CE, negative = BCE)",
    "lat": "Latitude",
    "lon": "Longitude",
}

# Data type requested for each coordinate variable by the data request.
GIAMIP_COORD_DTYPES = {
    "time": np.float64,
    "year": np.int32,
    "lat": np.float64,
    "lon": np.float64,
    "degree": np.float64,
    "order": np.float64,
}


def _same_text(actual: str, expected: str) -> bool:
    """Compare descriptive text ignoring case and runs of whitespace.

    long_name values are free text: the data request itself carries trailing
    spaces on several of them, so an exact string comparison would flag files
    that are in fact conformant. Wording differences are still reported.
    """
    return " ".join(str(actual).split()).casefold() == " ".join(expected.split()).casefold()


def _check_coord_dtype(log_file, ds: xr.Dataset, coord: str, kind: str) -> int:
    """Check one coordinate variable against the requested data type."""
    expected_dtype = GIAMIP_COORD_DTYPES.get(kind)
    if expected_dtype is None or ds[coord].dtype == expected_dtype:
        return 0

    log_file.write(
        f" - ERROR (coordinate): coordinate '{coord}' dtype is {ds[coord].dtype},"
        f" expected {np.dtype(expected_dtype).name}.\n"
    )
    return 1


def _check_coord_names(log_file, ds: xr.Dataset, coord: str, kind: str) -> int:
    """Check 'long_name', 'standard_name' and dtype on one coordinate variable.

    The data request asks for units and long_name on every coordinate variable,
    plus standard_name wherever a CF one applies, and gives a data type for each.
    *kind* selects the expected values from GIAMIP_COORD_LONG_NAMES /
    GIAMIP_COORD_STANDARD_NAMES / GIAMIP_COORD_DTYPES ("time", "year", "lat",
    "lon", ...), so a time axis named "t" is still checked against the values
    requested for time.
    """
    errors = _check_coord_dtype(log_file, ds, coord, kind)

    actual_ln = ds[coord].attrs.get("long_name")
    expected_ln = GIAMIP_COORD_LONG_NAMES.get(kind)
    if actual_ln is None:
        log_file.write(
            f" - ERROR (coordinate): coordinate '{coord}' missing 'long_name'.\n"
        )
        errors += 1
    elif expected_ln is not None and not _same_text(actual_ln, expected_ln):
        log_file.write(
            f" - ERROR (coordinate): coordinate '{coord}' long_name '{actual_ln}',"
            f" expected '{expected_ln}'.\n"
        )
        errors += 1

    expected_sn = GIAMIP_COORD_STANDARD_NAMES.get(kind)
    if expected_sn is not None:
        actual_sn = ds[coord].attrs.get("standard_name")
        if actual_sn is None:
            log_file.write(
                f" - ERROR (coordinate): coordinate '{coord}' missing 'standard_name'"
                f" (expected '{expected_sn}').\n"
            )
            errors += 1
        elif actual_sn != expected_sn:
            log_file.write(
                f" - ERROR (coordinate): coordinate '{coord}' standard_name"
                f" '{actual_sn}', expected '{expected_sn}'.\n"
            )
            errors += 1

    return errors


def _check_attributes(
    log_file,
    ds: xr.Dataset,
    var_name: str,
    var_meta: dict,
    is_spatial: bool,
    has_time: bool,
) -> tuple[int, int, int]:
    """Check attributes, split into (metadata, coordinate, variable) error counts.

    Returns a 3-tuple of error counts:
      - metadata:    global/file-level attributes (group, model, contributors).
      - coordinate:  attributes on the lat/lon/time coordinate variables.
      - variable:    attributes on the main data variable (units, long_name,
                     standard_name) and its dtype.
    """
    log_file.write("ATTRIBUTE Tests \n")

    expected_units = var_meta.get("units", "")

    # --- Metadata (global) attributes ---
    meta_errors = 0
    for attr in ("group", "model", "contributors_name", "contributors_email"):
        if attr not in ds.attrs:
            log_file.write(f" - ERROR (metadata): global attribute '{attr}' is missing.\n")
            meta_errors += 1

    if meta_errors == 0:
        log_file.write(" - Metadata attributes: passed\n")

    # --- Coordinate attributes ---
    coord_errors = 0
    if has_time:
        time_coord = _find_time_coord(ds)
        if time_coord is None:
            log_file.write(
                " - ERROR (coordinate): time coordinate variable not found.\n"
            )
            coord_errors += 1
        else:
            # units and calendar can sit in either attrs or encoding, depending on
            # how the file was written (datasets are opened with decode_times=False).
            combined = {**ds[time_coord].encoding, **ds[time_coord].attrs}
            u = combined.get("units")
            if u is None:
                log_file.write(
                    " - ERROR (coordinate): time coordinate missing 'units'.\n"
                )
                coord_errors += 1
            elif not _TIME_UNITS_RE.match(str(u)):
                log_file.write(
                    f" - ERROR (coordinate): time units '{u}',"
                    f" expected '{GIAMIP_TIME_UNITS}'.\n"
                )
                coord_errors += 1

            cal = combined.get("calendar")
            if cal is None:
                log_file.write(
                    " - ERROR (coordinate): time coordinate missing 'calendar'"
                    f" (expected '{GIAMIP_TIME_CALENDAR}').\n"
                )
                coord_errors += 1
            elif str(cal) != GIAMIP_TIME_CALENDAR:
                log_file.write(
                    f" - ERROR (coordinate): time calendar '{cal}',"
                    f" expected '{GIAMIP_TIME_CALENDAR}'.\n"
                )
                coord_errors += 1

            coord_errors += _check_coord_names(log_file, ds, time_coord, "time")

        # The companion 'year' variable. Its absence is reported by the time
        # checks, which need it; here we check it once it is there.
        if GIAMIP_YEAR_VAR in ds.variables:
            year_units = ds[GIAMIP_YEAR_VAR].attrs.get("units")
            if year_units != "year":
                log_file.write(
                    f" - ERROR (coordinate): '{GIAMIP_YEAR_VAR}' units"
                    f" '{year_units}', expected 'year'.\n"
                )
                coord_errors += 1
            coord_errors += _check_coord_names(
                log_file, ds, GIAMIP_YEAR_VAR, GIAMIP_YEAR_VAR
            )

    if is_spatial:
        for coord in ("lat", "lon"):
            if coord in ds.coords:
                if "units" not in ds[coord].attrs:
                    log_file.write(
                        f" - ERROR (coordinate): coordinate '{coord}' missing 'units'.\n"
                    )
                    coord_errors += 1
                coord_errors += _check_coord_names(log_file, ds, coord, coord)
            else:
                log_file.write(f" - ERROR (coordinate): coordinate '{coord}' not found.\n")
                coord_errors += 1

    # degree/order carry a requested data type too; their presence and the rest of
    # their attributes are not yet checked (see TODO.md).
    for coord in ("degree", "order"):
        if coord in ds.coords:
            coord_errors += _check_coord_dtype(log_file, ds, coord, coord)

    if coord_errors == 0:
        log_file.write(" - Coordinate attributes: passed\n")

    # --- Variable attributes ---
    var_errors = 0
    if var_name in ds:
        # Units
        actual_units = ds[var_name].attrs.get("units")
        if actual_units is None:
            if expected_units == "1":
                log_file.write(" - Units: absent (dimensionless '1' implied): passed\n")
            else:
                log_file.write(
                    f" - ERROR (variable): variable '{var_name}' has no 'units' attribute.\n"
                )
                var_errors += 1
        elif actual_units == expected_units:
            log_file.write(f" - Units '{actual_units}': passed\n")
        else:
            log_file.write(
                f" - ERROR (variable): units '{actual_units}', expected '{expected_units}'.\n"
            )
            var_errors += 1

        # long_name: must be present and match the wording of the data request.
        actual_ln = ds[var_name].attrs.get("long_name")
        expected_ln = var_meta.get("long_name")
        if actual_ln is None:
            log_file.write(
                f" - ERROR (variable): variable '{var_name}' missing 'long_name'.\n"
            )
            var_errors += 1
        elif expected_ln is not None and not _same_text(actual_ln, expected_ln):
            log_file.write(
                f" - ERROR (variable): variable '{var_name}' long_name '{actual_ln}'"
                f" does not match the data request '{expected_ln}'.\n"
            )
            var_errors += 1

        # standard_name (only checked when the data request specifies one)
        expected_sn = var_meta.get("standard_name")
        if expected_sn is not None:
            actual_sn = ds[var_name].attrs.get("standard_name")
            if actual_sn is None:
                log_file.write(
                    f" - ERROR (variable): variable '{var_name}' missing 'standard_name'.\n"
                )
                var_errors += 1
            elif actual_sn != expected_sn:
                log_file.write(
                    f" - ERROR (variable): standard_name '{actual_sn}'"
                    f" does not match expected '{expected_sn}'.\n"
                )
                var_errors += 1

        # dtype: per-variable precision from the data request (float32 for the
        # gridded fields, float64 for the scalar time series and Clm/Slm).
        expected_dtype = var_meta.get("precision")
        if expected_dtype is not None and ds[var_name].dtype != expected_dtype:
            log_file.write(
                f" - ERROR (variable): variable '{var_name}' dtype is {ds[var_name].dtype},"
                f" expected {np.dtype(expected_dtype).name}.\n"
            )
            var_errors += 1

        if var_errors == 0:
            log_file.write(" - Variable attributes: passed\n")
    else:
        # Reported as an error by _check_numerical; here the point is not to claim
        # that attribute checks passed when there was nothing to check them on.
        log_file.write(
            f" - Variable attributes: skipped (variable '{var_name}' is not in the file).\n"
        )

    return meta_errors, coord_errors, var_errors


def _run_variable_checks(
    log_file,
    ds: xr.Dataset,
    file_name: str,
    var_name: str,
    experiment_name: str,
    experiments: list,
    report_naming_issues: list,
    forcing=None,
) -> tuple[int, int, int, int, int, int, int]:
    naming_errors = 0
    num_errors = 0
    spatial_errors = 0
    dim_errors = 0
    time_errors = 0
    meta_attr_errors = 0
    coord_attr_errors = 0
    var_attr_errors = 0

    log_file.write(" \n")
    log_file.write(f"Experiment: {experiment_name} - File: {file_name}\n")
    log_file.write(" \n")

    naming_errors += _check_naming(
        log_file, file_name, var_name, experiment_name, experiments, report_naming_issues
    )
    if naming_errors:
        return (naming_errors, num_errors, spatial_errors, dim_errors, time_errors,
                meta_attr_errors, coord_attr_errors, var_attr_errors)

    var_meta = GIAMIP_VAR_META.get(var_name, {})
    dim_type = var_meta.get("dim", "t")
    is_spatial = dim_type == "lat_lon_t"
    has_time = dim_type in ("lat_lon_t", "t")

    log_file.write(f"** Tested Variable: {var_name}\n \n")

    forcing_years = forcing.get("years") if forcing else None
    forcing_lat = forcing.get("lat") if forcing else None
    forcing_lon = forcing.get("lon") if forcing else None

    num_errors += _check_numerical(log_file, ds, var_name)

    if is_spatial:
        spatial_errors += _check_spatial(
            log_file, ds, forcing_lat=forcing_lat, forcing_lon=forcing_lon
        )

    dim_errors += _check_dimensions(log_file, ds, var_name, dim_type)

    if has_time:
        time_errors += _check_time(
            log_file, ds,
            var_meta.get("output_interval", "forcing"),
            forcing_years=forcing_years,
        )

    m, c, v = _check_attributes(log_file, ds, var_name, var_meta, is_spatial, has_time)
    meta_attr_errors += m
    coord_attr_errors += c
    var_attr_errors += v

    return (naming_errors, num_errors, spatial_errors, dim_errors, time_errors,
            meta_attr_errors, coord_attr_errors, var_attr_errors)


# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------

def _print_experiment_summary(experiment_name: str, exp_errors: int) -> None:
    print(f"{experiment_name}: compliance check processed.")
    if exp_errors > 0:
        print(f"Found {exp_errors} error(s). Check compliance_checker_log.txt for details.")
    else:
        print("Successfully verified with no errors")
    print()


def _print_total_summary(source_path: str, total_errors: int) -> None:
    print("-------------------------------------------------------------------------")
    print(f"{source_path}: compliance check processed.")
    if total_errors > 0:
        print(f"Found a total of {total_errors} error(s). Check compliance_checker_log.txt.")
    else:
        print("Successfully verified with no errors")
    print("-------------------------------------------------------------------------")


def _strictly_increasing(values) -> bool:
    return all(x < y for x, y in zip(values, values[1:]))


def _write_log_header(log_file, commit_num: str, source_path: str, today: datetime.date) -> None:
    log_file.write(
        "************************************************************************************\n"
    )
    log_file.write(
        "*************     GIA Model Simulations - GIAMIP Compliance Checker    *************\n"
    )
    log_file.write(
        "************************************************************************************\n"
    )
    log_file.write(f"Commit Number: {commit_num} \n")
    log_file.write("verification criteria: GIAMIP data request\n")
    log_file.write(f"date: {today.strftime('%Y/%m/%d')}\n")
    log_file.write("source: https://github.com/JanJereczek/GIA_SimulationChecker \n")
    log_file.write(" \n")
    log_file.write(
        "------------------------------------------------------------------------------------\n"
    )
    log_file.write(f"Verified directory: {source_path} \n")
    log_file.write(
        "------------------------------------------------------------------------------------\n"
    )
    log_file.write(" \n")
    log_file.write(" \n")
    log_file.write(" \n")
    log_file.write(" \n")
    log_file.write(
        "====================================================================================\n"
    )
    log_file.write(
        "================                DETAILED RESULTS                    ================\n"
    )
    log_file.write(
        "====================================================================================\n"
    )
    log_file.write("Hint: Use Ctrl+F to look for specific problems. \n")
    log_file.write(" \n")


def _insert_synthesis(
    source_path: str,
    exp_counter: int,
    file_counter: int,
    total_errors: int,
    total_file_errors: int,
    total_naming_errors: int,
    total_num_errors: int,
    total_spatial_errors: int,
    total_dim_errors: int,
    total_time_errors: int,
    total_meta_attr_errors: int,
    total_coord_attr_errors: int,
    total_var_attr_errors: int,
    report_naming_issues: list,
) -> None:
    with open(os.path.join(source_path, "compliance_checker_log.txt"), "r") as f:
        contents = f.readlines()

    # Insert the synthesis into the blank region just after the "Verified directory"
    # header line, rather than relying on a hard-coded line number.
    header_idx = next(
        (i for i, line in enumerate(contents) if line.startswith("Verified directory")), 9
    )
    iline = header_idx + 2
    contents.insert(iline, f"{exp_counter} experiments checked.\n"); iline += 1
    contents.insert(iline, f"{file_counter} files checked.\n"); iline += 2
    contents.insert(iline, f"{total_errors} error(s) detected.\n"); iline += 1
    contents.insert(iline, f"  - Missing mandatory variable : {total_file_errors} error(s)\n"); iline += 1
    contents.insert(iline, f"  - Naming Tests               : {total_naming_errors} error(s)\n"); iline += 1
    contents.insert(iline, f"  - Numerical Tests            : {total_num_errors} error(s)\n"); iline += 1
    contents.insert(iline, f"  - Spatial Tests              : {total_spatial_errors} error(s)\n"); iline += 1
    contents.insert(iline, f"  - Dimension Tests            : {total_dim_errors} error(s)\n"); iline += 1
    contents.insert(iline, f"  - Time Tests                 : {total_time_errors} error(s)\n"); iline += 1
    contents.insert(iline, f"  - Metadata Attribute Tests   : {total_meta_attr_errors} error(s)\n"); iline += 1
    contents.insert(iline, f"  - Coordinate Attribute Tests : {total_coord_attr_errors} error(s)\n"); iline += 1
    contents.insert(iline, f"  - Variable Attribute Tests   : {total_var_attr_errors} error(s)\n"); iline += 2
    contents.insert(iline, "0 warning(s) detected.\n"); iline += 2

    if total_naming_errors > 0:
        contents.insert(iline, "Naming tests errors report: \n"); iline += 1
        for j, issue in enumerate(report_naming_issues):
            contents.insert(iline + j, f"  - {issue}\n")
        contents.insert(iline + len(report_naming_issues), "\n")

    with open(os.path.join(source_path, "compliance_checker_log.txt"), "w") as f:
        f.writelines(contents)


if __name__ == "__main__":
    main()
