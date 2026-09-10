"""
Pytest suite for the GIAMIP compliance checker.

Baseline: generate all mandatory variables for Exp01 with two time steps
(years 1 and 1001) on the full 257×513 global grid, then verify the checker
reports zero errors.  Each subsequent test mutates one aspect of the baseline
copy and verifies the checker catches that specific error.
"""

import importlib.util
import shutil
import sys
from pathlib import Path

import netCDF4
import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import giamip_compliance_checker as checker


def _load_generator():
    spec = importlib.util.spec_from_file_location(
        "generate_giamip_test_files",
        REPO_ROOT / "generate" / "generate_giamip_test_files.py",
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


generator = _load_generator()

EXP_ID = "Exp01"
GROUP = "TESTGROUP"
MODEL = "TESTMODEL"
START_YEAR = 1
END_YEAR = 1001


# ---------------------------------------------------------------------------
# Session-scoped baseline: generate all mandatory files once
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def baseline_dir(tmp_path_factory):
    root = tmp_path_factory.mktemp("giamip_baseline")
    generator.create_all_mandatory_files(
        output_dir=root,
        experiment_id=EXP_ID,
        group=GROUP,
        model=MODEL,
        start_year=START_YEAR,
        end_year=END_YEAR,
        n_steps=2,
    )
    summary = checker.run_checker(
        source_path=str(root),
        forcing_path=None,
        commit_num="tests",
    )
    assert summary["total_errors"] == 0, (
        "Baseline synthetic files should pass with zero errors, but found "
        f"{summary['total_errors']}.\n{summary['log_text']}"
    )
    return root


# Per-test copy of the baseline directory
@pytest.fixture
def case_dir(tmp_path, baseline_dir):
    case = tmp_path / "CORE"
    shutil.copytree(baseline_dir, case)
    return case


def run_checker(case_dir: Path) -> dict:
    return checker.run_checker(
        source_path=str(case_dir),
        forcing_path=None,
        commit_num="tests",
    )


# Encoding for the 'time'/'year' coordinate pair, as the data request asks for it.
TIME_ENCODING = {
    "time": {"dtype": "f8", "_FillValue": None},
    checker.GIAMIP_YEAR_VAR: {"dtype": "i4", "_FillValue": None},
}


def first_file(case_dir: Path) -> Path:
    """Return the first .nc file in the case directory."""
    return sorted(case_dir.glob("*.nc"))[0]


def file_section(log_text: str, file_name: str) -> str:
    """Return the part of the log covering a single file.

    The log concatenates every file in the directory, so a bare `in log_text`
    assertion can be satisfied by a different file's section. Each file's block
    starts with its "Experiment: ... - File: <name>" banner.
    """
    marker = f"- File: {file_name}\n"
    start = log_text.index(marker) + len(marker)
    nxt = log_text.find("- File: ", start)
    return log_text[start:] if nxt == -1 else log_text[start:nxt]


def file_for(case_dir: Path, var_name: str) -> Path:
    """Return the file for the given variable."""
    matches = list(case_dir.glob(f"{var_name}_*.nc"))
    assert matches, f"No file found for variable '{var_name}'"
    return matches[0]


# ---------------------------------------------------------------------------
# Helpers for mutating files
# ---------------------------------------------------------------------------

def set_global_attr(path: Path, attr: str, value) -> None:
    with netCDF4.Dataset(path, "a") as nc:
        setattr(nc, attr, value)


def del_global_attr(path: Path, attr: str) -> None:
    with netCDF4.Dataset(path, "a") as nc:
        nc.delncattr(attr)


def set_var_attr(path: Path, var: str, attr: str, value) -> None:
    with netCDF4.Dataset(path, "a") as nc:
        setattr(nc.variables[var], attr, value)


def del_var_attr(path: Path, var: str, attr: str) -> None:
    with netCDF4.Dataset(path, "a") as nc:
        nc.variables[var].delncattr(attr)


def inject_nan(path: Path, var: str) -> None:
    """Set one value in *var* to NaN."""
    with netCDF4.Dataset(path, "a") as nc:
        v = nc.variables[var]
        data = v[:]
        flat = data.flatten()
        flat[0] = np.nan
        v[:] = flat.reshape(data.shape)


def set_time_years(path: Path, years: list[int]) -> None:
    """Replace the time axis with the given years, keeping 'time' and 'year' consistent."""
    years = np.asarray(years, dtype=np.int32)
    with netCDF4.Dataset(path, "a") as nc:
        nc.variables["time"][:] = years.astype(np.float64) * checker.GIAMIP_DAYS_PER_YEAR
        nc.variables[checker.GIAMIP_YEAR_VAR][:] = years


def write_forcing(path: Path, lat, lon, years) -> Path:
    """Write a minimal forcing file with a 'year' variable and a lat/lon grid."""
    import xarray as xr

    ds = xr.Dataset(
        {"year": ("time", np.asarray(years, dtype=np.int32))},
        coords={
            "time": ("time", np.arange(len(years), dtype=np.float32)),
            "lat": ("lat", np.asarray(lat, dtype=np.float64)),
            "lon": ("lon", np.asarray(lon, dtype=np.float64)),
        },
    )
    ds.to_netcdf(path)
    return path


# ---------------------------------------------------------------------------
# Tests: baseline
# ---------------------------------------------------------------------------

def test_baseline_passes(case_dir):
    summary = run_checker(case_dir)
    assert summary["total_errors"] == 0
    assert "No errors. Good job !" in summary["log_text"]


# ---------------------------------------------------------------------------
# Tests: naming errors
# ---------------------------------------------------------------------------

def test_wrong_number_of_filename_fields(case_dir):
    f = first_file(case_dir)
    bad_name = f.stem + "_EXTRA.nc"
    f.rename(f.parent / bad_name)

    summary = run_checker(case_dir)

    # File cannot be parsed → grouped under _unknown → experiment not recognised
    assert summary["total_naming_errors"] >= 1
    assert summary["total_errors"] >= 1
    assert "_unknown" in summary["log_text"] or "not in" in summary["log_text"]


def test_unknown_experiment_in_filename(case_dir):
    # Rename every file so the experiment field is "UNKNOWN"
    for nc_file in list(case_dir.glob("*.nc")):
        parts = nc_file.stem.split("_")
        parts[checker.GIAMIP_FILENAME_EXP_IDX] = "UNKNOWN"
        nc_file.rename(nc_file.parent / ("_".join(parts) + ".nc"))

    summary = run_checker(case_dir)

    assert summary["total_errors"] >= 1
    assert "not in" in summary["log_text"]


def test_missing_mandatory_variable(case_dir):
    file_for(case_dir, "maf").unlink()

    summary = run_checker(case_dir)

    assert summary["total_file_errors"] >= 1
    assert "mandatory variable(s) are missing" in summary["log_text"]


# ---------------------------------------------------------------------------
# Tests: filename parsing, incl. optional ensemble-member identifier (_mYY)
# ---------------------------------------------------------------------------

def test_parse_filename_without_member():
    assert checker._parse_giamip_filename("delta_g_Exp03_JPL_ISSM-SLC.nc") == (
        "delta_g", "Exp03", "JPL", "ISSM-SLC"
    )


def test_parse_filename_with_member():
    # Member is folded into the returned experiment id.
    assert checker._parse_giamip_filename("delta_g_Exp10_m01_JPL_ISSM-SLC.nc") == (
        "delta_g", "Exp10_m01", "JPL", "ISSM-SLC"
    )
    # Works for multi-token variable names too.
    assert checker._parse_giamip_filename(
        "ocean_area_fraction_Exp08_m02_GRP_MOD.nc"
    ) == ("ocean_area_fraction", "Exp08_m02", "GRP", "MOD")


def test_parse_filename_members_are_distinct():
    m01 = checker._parse_giamip_filename("maf_Exp10_m01_JPL_ISSM-SLC.nc")
    m02 = checker._parse_giamip_filename("maf_Exp10_m02_JPL_ISSM-SLC.nc")
    assert m01[1] == "Exp10_m01"
    assert m02[1] == "Exp10_m02"
    assert m01[1] != m02[1]


def test_parse_filename_extra_field_is_rejected():
    # A 4th trailing field that is not a valid member (mYY) is malformed.
    assert checker._parse_giamip_filename("maf_Exp10_NOTAMEMBER_JPL_MOD.nc") is None


def test_ensemble_members_treated_as_separate_experiments(tmp_path):
    """m01 and m02 are independent experiments, each requiring all mandatory vars."""
    root = tmp_path / "members"
    for member in ("Exp01_m01", "Exp01_m02"):
        generator.create_all_mandatory_files(
            output_dir=root, experiment_id=member, group=GROUP, model=MODEL,
            start_year=START_YEAR, end_year=END_YEAR, n_steps=2,
        )
    # Drop one mandatory variable from m02 only.
    next(root.glob("maf_Exp01_m02_*.nc")).unlink()

    summary = checker.run_checker(
        source_path=str(root), forcing_path=None,
        commit_num="tests",
    )

    # Only m02 should report a missing mandatory variable; m01 is complete.
    assert summary["total_file_errors"] >= 1
    assert "In experiment Exp01_m02" in summary["log_text"]
    assert "Exp01_m01: all mandatory variables present" in summary["log_text"]


# ---------------------------------------------------------------------------
# Tests: model- and full-level orchestration (loop over experiments / models)
# ---------------------------------------------------------------------------

def test_forcing_mapping():
    # Exp05 and Exp08 share the PaleoMIST forcing.
    assert checker.EXPERIMENT_FORCING["Exp05"] == "iceHistory-PaleoMIST_1a.nc"
    assert checker.EXPERIMENT_FORCING["Exp05"] == checker.EXPERIMENT_FORCING["Exp08"]
    # All 12 experiments are mapped.
    assert set(checker.EXPERIMENT_FORCING) == set(checker.GIAMIP_EXPERIMENTS)
    # Member ids resolve via the base experiment.
    p = checker._forcing_path_for_experiment("Exp05_m02", "/some/input")
    assert p.endswith("iceHistory-PaleoMIST_1a.nc")
    assert "some/input" in p
    # Unknown experiment -> no mapping.
    assert checker._forcing_path_for_experiment("NotAnExp", "/x") is None


def _make_experiment_dir(parent, experiment_id):
    generator.create_all_mandatory_files(
        output_dir=parent / experiment_id, experiment_id=experiment_id,
        group=GROUP, model=MODEL, start_year=START_YEAR, end_year=END_YEAR, n_steps=2,
    )


def test_run_model_checker_loops_over_experiments(tmp_path):
    model = tmp_path / "SomeModel"
    _make_experiment_dir(model, "Exp01")
    _make_experiment_dir(model, "Exp05")
    forcing_dir = tmp_path / "input"  # intentionally empty (forcing not delivered yet)
    forcing_dir.mkdir()

    results = checker.run_model_checker(
        source_dir=str(model), forcing_dir=str(forcing_dir),
        commit_num="tests",
    )

    assert set(results) == {"Exp01", "Exp05"}
    # No forcing files present -> only intrinsic checks; generated files comply.
    assert results["Exp01"]["total_errors"] == 0
    assert results["Exp05"]["total_errors"] == 0


def test_run_model_checker_skips_non_experiment_dirs(tmp_path):
    model = tmp_path / "SomeModel"
    _make_experiment_dir(model, "Exp05")
    (model / "notes").mkdir()  # should be ignored
    forcing_dir = tmp_path / "input"
    forcing_dir.mkdir()

    results = checker.run_model_checker(
        source_dir=str(model), forcing_dir=str(forcing_dir),
        commit_num="tests",
    )

    assert set(results) == {"Exp05"}


def test_run_full_checker_loops_over_models(tmp_path):
    out = tmp_path / "output"
    for model_name in ("ModelA", "ModelB"):
        _make_experiment_dir(out / model_name, "Exp05")
    forcing_dir = tmp_path / "input"
    forcing_dir.mkdir()

    results = checker.run_full_checker(
        output_dir=str(out), forcing_dir=str(forcing_dir),
        commit_num="tests",
    )

    assert set(results) == {"ModelA", "ModelB"}
    assert results["ModelA"]["Exp05"]["total_errors"] == 0
    assert results["ModelB"]["Exp05"]["total_errors"] == 0


# ---------------------------------------------------------------------------
# Tests: numerical errors
# ---------------------------------------------------------------------------

def test_nan_values_detected(case_dir):
    path = file_for(case_dir, "maf")
    inject_nan(path, "maf")

    summary = run_checker(case_dir)

    assert summary["total_num_errors"] >= 1
    assert "NaN" in summary["log_text"] or "missing" in summary["log_text"]


def test_wrong_units(case_dir):
    path = file_for(case_dir, "delta_bed")
    set_var_attr(path, "delta_bed", "units", "km")

    summary = run_checker(case_dir)

    assert summary["total_var_attr_errors"] >= 1
    assert "units" in summary["log_text"]


def test_mask_out_of_range(case_dir):
    """Inject values > 1 into a mask variable."""
    path = file_for(case_dir, "ocean_area_fraction")
    with netCDF4.Dataset(path, "a") as nc:
        nc.variables["ocean_area_fraction"][0, 0, 0] = 1.5

    summary = run_checker(case_dir)

    assert summary["total_num_errors"] >= 1
    assert "outside its plausible range" in summary["log_text"]


def test_out_of_bounds_detected(case_dir):
    """A scalar value outside the variable's plausible range is flagged."""
    path = file_for(case_dir, "grd_ice_mass")
    with netCDF4.Dataset(path, "a") as nc:
        nc.variables["grd_ice_mass"][0] = 1e25  # far above the 1e20 upper bound

    summary = run_checker(case_dir)

    assert summary["total_num_errors"] >= 1
    assert "outside its plausible range" in summary["log_text"]


# ---------------------------------------------------------------------------
# Tests: dimension ordering
# ---------------------------------------------------------------------------

def test_baseline_has_no_dimension_errors(case_dir):
    summary = run_checker(case_dir)
    assert summary["total_dim_errors"] == 0
    assert "Dimension ordering" in summary["log_text"]


def test_wrong_dimension_order_detected(case_dir):
    """Transpose a field to (lat, lon, time) — wrong order — and expect an error."""
    import xarray as xr

    path = file_for(case_dir, "delta_bed")
    with xr.open_dataset(path, decode_times=False) as ds:
        ds = ds.load()
    ds["delta_bed"] = ds["delta_bed"].transpose("lat", "lon", "time")
    path.unlink()
    ds.to_netcdf(path)

    summary = run_checker(case_dir)

    assert summary["total_dim_errors"] >= 1
    assert "dimension ordering" in summary["log_text"].lower()
    # (lat, lon, time) is not the reverse of (time, lat, lon), so no column-major hint
    assert "HINT" not in summary["log_text"]


def test_reversed_dimension_order_hints_at_column_major(case_dir):
    """(lon, lat, time) is what Fortran/MATLAB produce from a literal transcription."""
    import xarray as xr

    path = file_for(case_dir, "delta_bed")
    with xr.open_dataset(path, decode_times=False) as ds:
        ds = ds.load()
    ds["delta_bed"] = ds["delta_bed"].transpose("lon", "lat", "time")
    path.unlink()
    ds.to_netcdf(path)

    summary = run_checker(case_dir)

    assert summary["total_dim_errors"] >= 1
    assert "HINT: the dimensions are in exactly the reverse" in summary["log_text"]
    assert "ncdump -h" in summary["log_text"]


def test_reversed_degree_order_hints_at_column_major(case_dir):
    """Same trap for Clm/Slm, where both axes are 97 long so sizes give nothing away."""
    import xarray as xr

    generator.create_giamip_file(
        output_dir=case_dir, variable_name="Clm", experiment_id=EXP_ID,
        group=GROUP, model=MODEL,
    )
    path = file_for(case_dir, "Clm")
    with xr.open_dataset(path, decode_times=False) as ds:
        ds = ds.load()
    ds["Clm"] = ds["Clm"].transpose("order", "degree")
    path.unlink()
    ds.to_netcdf(path)

    summary = run_checker(case_dir)

    assert summary["total_dim_errors"] >= 1
    assert "HINT: the dimensions are in exactly the reverse" in summary["log_text"]


# ---------------------------------------------------------------------------
# Tests: spatial errors
# ---------------------------------------------------------------------------

def test_wrong_grid_size(case_dir, tmp_path):
    """Replace a 3D file with one that has a 10×20 grid instead of 257×513."""
    path = file_for(case_dir, "delta_bed")
    path.unlink()

    # Generate a small-grid replacement
    small_lat = np.linspace(-81.0, 81.0, 10, dtype=np.float64)
    small_lon = np.linspace(0.0, 360.0, 20, endpoint=False, dtype=np.float64)
    time_years = np.array([START_YEAR, END_YEAR], dtype=np.int32)
    data = np.random.uniform(0.0, 100.0, (2, 10, 20)).astype(np.float32)

    ds = _make_spatial_dataset("delta_bed", data, time_years, small_lat, small_lon)
    ds.to_netcdf(path, unlimited_dims=("time",),
                 encoding={"delta_bed": {"dtype": "f4", "_FillValue": None},
                            **TIME_ENCODING})

    summary = run_checker(case_dir)

    assert summary["total_spatial_errors"] >= 1
    assert "nodes, expected" in summary["log_text"]


def test_wrong_lat_range(case_dir):
    """Replace a 3D file with one where lat spans [0, 180] instead of [-90, 90]."""
    path = file_for(case_dir, "delta_g")
    path.unlink()

    bad_lat = np.linspace(0.0, 180.0, 257, dtype=np.float64)
    lon = np.linspace(0.0, 360.0, 513, endpoint=False, dtype=np.float64)
    time_years = np.array([START_YEAR, END_YEAR], dtype=np.int32)
    data = np.random.uniform(-10.0, 10.0, (2, 257, 513)).astype(np.float32)

    ds = _make_spatial_dataset("delta_g", data, time_years, bad_lat, lon)
    ds.to_netcdf(path, unlimited_dims=("time",),
                 encoding={"delta_g": {"dtype": "f4", "_FillValue": None},
                            **TIME_ENCODING})

    summary = run_checker(case_dir)

    assert summary["total_spatial_errors"] >= 1
    assert "north edge" in summary["log_text"]


def test_forcing_grid_match_passes(case_dir, tmp_path):
    """A forcing file whose grid matches the data passes the grid comparison."""
    lat, lon = generator._make_grid()
    forcing = write_forcing(
        tmp_path / "forcing.nc", lat, lon, [START_YEAR, END_YEAR]
    )

    summary = checker.run_checker(
        source_path=str(case_dir),
        forcing_path=str(forcing),
        commit_num="tests",
    )

    assert summary["total_spatial_errors"] == 0
    assert "grid matches forcing: passed" in summary["log_text"]


def test_forcing_grid_mismatch_detected(case_dir, tmp_path):
    """A forcing file with a shifted latitude grid is flagged as a spatial error."""
    lat, lon = generator._make_grid()
    shifted_lat = lat + 0.05  # well beyond the relative tolerance
    forcing = write_forcing(
        tmp_path / "forcing.nc", shifted_lat, lon, [START_YEAR, END_YEAR]
    )

    summary = checker.run_checker(
        source_path=str(case_dir),
        forcing_path=str(forcing),
        commit_num="tests",
    )

    assert summary["total_spatial_errors"] >= 1
    assert "does not match the forcing grid" in summary["log_text"]


# ---------------------------------------------------------------------------
# Tests: time errors
# ---------------------------------------------------------------------------

def test_non_monotonic_time(case_dir):
    """Reverse the time axis so it is decreasing."""
    path = file_for(case_dir, "grd_ice_mass")
    set_time_years(path, [END_YEAR, START_YEAR])

    summary = run_checker(case_dir)

    assert summary["total_time_errors"] >= 1
    assert "monotonically" in summary["log_text"]


# ---------------------------------------------------------------------------
# Tests: attribute errors
# ---------------------------------------------------------------------------

def test_missing_global_attribute(case_dir):
    del_global_attr(first_file(case_dir), "contributors_email")

    summary = run_checker(case_dir)

    assert summary["total_meta_attr_errors"] >= 1
    assert "contributors_email" in summary["log_text"]


def test_wrong_time_units(case_dir):
    """The pre-v3 convention (calendar years in 'time') must be flagged."""
    set_var_attr(first_file(case_dir), "time", "units", "year")

    summary = run_checker(case_dir)

    assert summary["total_coord_attr_errors"] >= 1
    assert "expected 'days since 2000-01-01 00:00:00'" in summary["log_text"]


def test_short_time_units_accepted(case_dir):
    """CF treats 'days since 2000-01-01' as identical; it must not be flagged."""
    for nc_file in case_dir.glob("*.nc"):
        with netCDF4.Dataset(nc_file, "a") as nc:
            if "time" in nc.variables:
                nc.variables["time"].units = "days since 2000-01-01"

    summary = run_checker(case_dir)

    assert summary["total_errors"] == 0, summary["log_text"]


def test_wrong_time_calendar(case_dir):
    set_var_attr(first_file(case_dir), "time", "calendar", "standard")

    summary = run_checker(case_dir)

    assert summary["total_coord_attr_errors"] >= 1
    assert "expected '365_day'" in summary["log_text"]


def test_missing_year_variable(case_dir):
    path = file_for(case_dir, "maf")
    _drop_variable(path, checker.GIAMIP_YEAR_VAR)

    summary = run_checker(case_dir)

    assert summary["total_time_errors"] >= 1
    assert "'year' not found" in summary["log_text"]


def test_year_inconsistent_with_time(case_dir):
    """Years written into 'time' as well as 'year' — the old convention."""
    path = file_for(case_dir, "maf")
    with netCDF4.Dataset(path, "a") as nc:
        nc.variables["time"][:] = np.array([START_YEAR, END_YEAR], dtype=np.float64)

    summary = run_checker(case_dir)

    assert summary["total_time_errors"] >= 1
    assert "does not match 'time'/365" in summary["log_text"]


def test_wrong_coordinate_dtype(case_dir):
    """lat is requested as float64."""
    path = file_for(case_dir, "delta_bed")
    _rewrite_delta_bed(path, lat_dtype="f4")

    summary = run_checker(case_dir)

    assert summary["total_coord_attr_errors"] >= 1
    assert "'lat' dtype is float32, expected float64" in summary["log_text"]


def test_missing_long_name(case_dir):
    path = file_for(case_dir, "maf")
    del_var_attr(path, "maf", "long_name")

    summary = run_checker(case_dir)

    assert summary["total_var_attr_errors"] >= 1
    assert "long_name" in summary["log_text"]


def test_wrong_dtype(case_dir, tmp_path):
    """Replace the delta_bed file with one where the variable is float64."""
    path = file_for(case_dir, "delta_bed")
    _rewrite_delta_bed(path, dtype="f8")

    summary = run_checker(case_dir)

    assert summary["total_var_attr_errors"] >= 1
    assert "expected float32" in summary["log_text"]


def test_wrong_dtype_scalar(case_dir):
    """Scalar time series are requested as float64; float32 must be flagged."""
    path = file_for(case_dir, "maf")
    _rewrite_scalar(path, "maf", dtype="f4")

    summary = run_checker(case_dir)

    assert summary["total_var_attr_errors"] >= 1
    assert "expected float64" in summary["log_text"]


def test_wrong_long_name(case_dir):
    path = file_for(case_dir, "maf")
    set_var_attr(path, "maf", "long_name", "Mass above flotation")

    summary = run_checker(case_dir)

    assert summary["total_var_attr_errors"] >= 1
    assert "does not match the data request" in summary["log_text"]


def test_long_name_whitespace_and_case_tolerated(case_dir):
    """The data request itself carries stray spaces, so those must not fail."""
    path = file_for(case_dir, "maf")
    expected = checker.GIAMIP_VAR_META["maf"]["long_name"]
    set_var_attr(path, "maf", "long_name", f"  {expected.upper()}  ")

    summary = run_checker(case_dir)

    assert summary["total_errors"] == 0, summary["log_text"]


def test_wrong_coordinate_long_name(case_dir):
    set_var_attr(file_for(case_dir, "delta_bed"), "lat", "long_name", "y")

    summary = run_checker(case_dir)

    assert summary["total_coord_attr_errors"] >= 1
    assert "expected 'Latitude'" in summary["log_text"]


def test_missing_coordinate_long_name(case_dir):
    del_var_attr(file_for(case_dir, "delta_bed"), "lat", "long_name")

    summary = run_checker(case_dir)

    assert summary["total_coord_attr_errors"] >= 1
    assert "'lat' missing 'long_name'" in summary["log_text"]


def test_missing_coordinate_standard_name(case_dir):
    del_var_attr(file_for(case_dir, "delta_bed"), "lon", "standard_name")

    summary = run_checker(case_dir)

    assert summary["total_coord_attr_errors"] >= 1
    assert "'lon' missing 'standard_name'" in summary["log_text"]


def test_wrong_coordinate_standard_name(case_dir):
    set_var_attr(file_for(case_dir, "delta_bed"), "time", "standard_name", "t")

    summary = run_checker(case_dir)

    assert summary["total_coord_attr_errors"] >= 1
    assert "expected 'time'" in summary["log_text"]


# ---------------------------------------------------------------------------
# Tests: declared fill values
# ---------------------------------------------------------------------------

def test_declared_fill_value(case_dir):
    """_FillValue is declared in the file and lands in xarray's encoding."""
    path = file_for(case_dir, "delta_bed")
    _rewrite_delta_bed(path, fill_value=9.96921e36)

    summary = run_checker(case_dir)

    assert summary["total_num_errors"] >= 1
    assert "_FillValue" in summary["log_text"]


def test_declared_missing_value(case_dir):
    set_var_attr(file_for(case_dir, "maf"), "maf", "missing_value", -9999.0)

    summary = run_checker(case_dir)

    assert summary["total_num_errors"] >= 1
    assert "missing_value" in summary["log_text"]


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _drop_variable(path: Path, name: str) -> None:
    """Rewrite a file without *name* (netCDF4 cannot delete variables in place)."""
    import xarray as xr

    ds = xr.open_dataset(path, decode_times=False).load()
    ds = ds.drop_vars(name)
    encoding = {k: v for k, v in TIME_ENCODING.items() if k in ds.variables}
    path.unlink()
    ds.to_netcdf(path, unlimited_dims=("time",), encoding=encoding)


def _rewrite_delta_bed(path: Path, dtype: str = "f4", fill_value=None,
                       lat_dtype: str = "f8") -> None:
    """Rewrite the delta_bed file with the given on-disk dtypes and _FillValue.

    All three are encoding-level properties, so the file has to be recreated
    rather than patched in place.
    """
    path.unlink()

    lat, lon = generator._make_grid()
    time_years = np.array([START_YEAR, END_YEAR], dtype=np.int32)
    data = np.random.uniform(-100.0, 100.0, (2, 257, 513)).astype(dtype)

    ds = _make_spatial_dataset("delta_bed", data, time_years, lat, lon)
    ds.to_netcdf(
        path, unlimited_dims=("time",),
        encoding={"delta_bed": {"dtype": dtype, "_FillValue": fill_value},
                  "lat": {"dtype": lat_dtype, "_FillValue": None},
                  **TIME_ENCODING},
    )


def _rewrite_scalar(path: Path, var_name: str, dtype: str) -> None:
    """Rewrite a scalar (time,) variable file with the given on-disk dtype."""
    import xarray as xr
    from giamip_compliance_checker import GIAMIP_VAR_META

    meta = GIAMIP_VAR_META[var_name]
    lo, hi = meta["bounds"]
    values = np.linspace(lo + 0.1 * (hi - lo), hi - 0.1 * (hi - lo), 2).astype(dtype)
    time_years = np.array([START_YEAR, END_YEAR], dtype=np.int32)

    attrs = {"long_name": meta["long_name"], "units": meta["units"]}
    if meta.get("standard_name"):
        attrs["standard_name"] = meta["standard_name"]

    ds = xr.Dataset(
        {var_name: (("time",), values, attrs)},
        coords=dict(generator._time_coords(time_years)),
    )
    ds.attrs.update({
        "group": GROUP, "model": MODEL,
        "contributors_name": "Synthetic Test",
        "contributors_email": "test@example.org",
    })

    path.unlink()
    ds.to_netcdf(
        path, unlimited_dims=("time",),
        encoding={var_name: {"dtype": dtype, "_FillValue": None},
                  **TIME_ENCODING},
    )


def _make_spatial_dataset(var_name, data, time_years, lat, lon):
    """Build a minimal xarray Dataset for a (time, lat, lon) variable."""
    from giamip_compliance_checker import GIAMIP_VAR_META
    import xarray as xr

    meta = GIAMIP_VAR_META[var_name]
    attrs = {"long_name": meta["long_name"], "units": meta["units"]}
    if meta.get("standard_name"):
        attrs["standard_name"] = meta["standard_name"]

    ds = xr.Dataset(
        {var_name: (("time", "lat", "lon"), data, attrs)},
        coords={
            **generator._time_coords(time_years),
            "lat": ("lat", lat, generator.COORD_ATTRS["lat"]),
            "lon": ("lon", lon, generator.COORD_ATTRS["lon"]),
        },
    )
    ds.attrs.update({
        "group": GROUP, "model": MODEL,
        "contributors_name": "Synthetic Test",
        "contributors_email": "test@example.org",
    })
    return ds


# ---------------------------------------------------------------------------
# Tests: a directory that cannot be checked must never report zero errors
# ---------------------------------------------------------------------------

def test_empty_directory_is_not_compliant(tmp_path):
    """An experiment that was never delivered must not pass."""
    empty = tmp_path / "Exp01"
    empty.mkdir()

    summary = checker.run_checker(
        source_path=str(empty), forcing_path=None,
        commit_num="tests",
    )

    assert summary["total_errors"] == len(checker.GIAMIP_MANDATORY_VARS)
    assert summary["total_file_errors"] == len(checker.GIAMIP_MANDATORY_VARS)
    assert "No .nc files found" in summary["log_text"]
    # The synthesis block is written for this case too.
    assert "error(s) detected." in summary["log_text"]


def test_missing_directory_is_not_compliant(tmp_path):
    summary = checker.run_checker(
        source_path=str(tmp_path / "does_not_exist"), forcing_path=None,
        commit_num="tests",
    )

    assert summary["total_errors"] == len(checker.GIAMIP_MANDATORY_VARS)


def test_model_rollup_flags_empty_experiment(tmp_path):
    """The roll-up must not show '0 error(s)' for an experiment with no files."""
    model = tmp_path / "SomeModel"
    _make_experiment_dir(model, "Exp05")
    (model / "Exp01").mkdir()  # delivered as an empty directory
    forcing_dir = tmp_path / "input"
    forcing_dir.mkdir()

    results = checker.run_model_checker(
        source_dir=str(model), forcing_dir=str(forcing_dir),
        commit_num="tests",
    )

    assert results["Exp05"]["total_errors"] == 0
    assert results["Exp01"]["total_errors"] > 0


# ---------------------------------------------------------------------------
# Tests: the prescribed grid is enforced without a forcing file
# ---------------------------------------------------------------------------

def test_canonical_grid_shape_and_extent():
    lat = checker._canonical_lat()
    lon = checker._canonical_lon()

    assert lat.shape == (checker.GIAMIP_GRID_NLAT,)
    assert lon.shape == (checker.GIAMIP_GRID_NLON,)
    # Gaussian latitudes never reach the poles, and 360 is not repeated as 0.
    assert -90.0 < lat[0] < -89.0
    assert 89.0 < lat[-1] < 90.0
    assert lon[0] == 0.0
    assert 359.0 < lon[-1] < 360.0
    # Latitude spacing varies (Gaussian); longitude spacing does not.
    assert not np.allclose(np.diff(lat), np.diff(lat)[0], atol=1e-6)
    assert np.allclose(np.diff(lon), np.diff(lon)[0])


def _replace_grid(case_dir: Path, var_name: str, lat, lon) -> None:
    """Rewrite *var_name* on the given lat/lon axes, leaving everything else alone."""
    import xarray as xr

    path = file_for(case_dir, var_name)
    with xr.open_dataset(path, decode_times=False) as ds:
        ds = ds.load()
    ds = ds.assign_coords(
        lat=("lat", np.asarray(lat, dtype=np.float64), ds["lat"].attrs),
        lon=("lon", np.asarray(lon, dtype=np.float64), ds["lon"].attrs),
    )
    path.unlink()
    ds.to_netcdf(
        path, unlimited_dims=("time",),
        encoding={var_name: {"dtype": "f4", "_FillValue": None}, **TIME_ENCODING},
    )


def test_regional_grid_detected_without_forcing(case_dir):
    """A 257x513 grid covering only a region must not pass, forcing file or not."""
    _replace_grid(
        case_dir, "delta_bed",
        lat=np.linspace(-30.0, 30.0, 257),
        lon=np.linspace(0.0, 40.0, 513),
    )

    summary = run_checker(case_dir)

    assert summary["total_spatial_errors"] >= 1
    assert "does not match the prescribed GIAMIP grid" in summary["log_text"]
    assert "looks like a regional grid" in summary["log_text"]


def test_equally_spaced_latitude_detected_without_forcing(case_dir):
    """A global but equally spaced latitude axis is not the Gaussian grid."""
    _, lon = generator._make_grid()
    _replace_grid(case_dir, "delta_bed", lat=np.linspace(-89.6, 89.6, 257), lon=lon)

    summary = run_checker(case_dir)

    assert summary["total_spatial_errors"] >= 1
    assert "the latitude axis is equally spaced" in summary["log_text"]
    # Longitude is untouched and must still pass.
    assert "Longitude matches the GIAMIP grid: passed" in summary["log_text"]


# ---------------------------------------------------------------------------
# Tests: a check that could not run says "skipped", never "passed"
# ---------------------------------------------------------------------------

def test_absent_variable_does_not_report_attributes_passed(case_dir):
    """A file whose data variable is missing must not claim its attributes passed."""
    import xarray as xr

    path = file_for(case_dir, "mean_delta_g")
    with xr.open_dataset(path, decode_times=False) as ds:
        ds = ds.load()
    ds = ds.rename({"mean_delta_g": "mdg"})
    path.unlink()
    ds.to_netcdf(path, unlimited_dims=("time",), encoding=TIME_ENCODING)

    summary = run_checker(case_dir)
    section = file_section(summary["log_text"], path.name)

    assert "variable 'mean_delta_g' not found in dataset" in section
    assert "Variable attributes: skipped (variable 'mean_delta_g' is not in the file)" in section
    assert "Dimension ordering: skipped (variable 'mean_delta_g' is not in the file)" in section
    assert "Variable attributes: passed" not in section


# ---------------------------------------------------------------------------
# Tests: a missing time axis must not be validated against a fabricated one
# ---------------------------------------------------------------------------

def test_missing_time_coordinate_detected(case_dir):
    """xarray invents an index for a bare dimension; the checker must not trust it."""
    import xarray as xr

    path = file_for(case_dir, "maf")
    with xr.open_dataset(path, decode_times=False) as ds:
        ds = ds.load()
    # Drop the coordinate variables but keep the 'time' dimension on 'maf'.
    ds = ds.drop_vars(["time", checker.GIAMIP_YEAR_VAR])
    path.unlink()
    ds.to_netcdf(path)

    summary = run_checker(case_dir)
    section = file_section(summary["log_text"], path.name)

    assert summary["total_time_errors"] >= 1
    assert "no 'time' coordinate variable" in section
    # The fabricated 0, 1, 2, ... index must not be reported as a valid time axis.
    assert "Time is monotonically increasing: passed" not in section
