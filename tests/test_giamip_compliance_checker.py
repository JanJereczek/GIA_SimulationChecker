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
        workdir=str(REPO_ROOT),
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
        workdir=str(REPO_ROOT),
        commit_num="tests",
    )


def first_file(case_dir: Path) -> Path:
    """Return the first .nc file in the case directory."""
    return sorted(case_dir.glob("*.nc"))[0]


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
    """Replace the time axis with the given calendar years (stored directly)."""
    with netCDF4.Dataset(path, "a") as nc:
        nc.variables["time"][:] = np.array(years, dtype=np.float32)


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
# Tests: spatial errors
# ---------------------------------------------------------------------------

def test_wrong_grid_size(case_dir, tmp_path):
    """Replace a 3D file with one that has a 10×20 grid instead of 257×513."""
    path = file_for(case_dir, "delta_bed")
    path.unlink()

    # Generate a small-grid replacement
    small_lat = np.linspace(-81.0, 81.0, 10, dtype=np.float32)
    small_lon = np.linspace(0.0, 360.0, 20, endpoint=False, dtype=np.float32)
    time_years = np.array([START_YEAR, END_YEAR], dtype=np.float32)
    data = np.random.uniform(0.0, 100.0, (2, 10, 20)).astype(np.float32)

    ds = _make_spatial_dataset("delta_bed", data, time_years, small_lat, small_lon)
    ds.to_netcdf(path, unlimited_dims=("time",),
                 encoding={"delta_bed": {"dtype": "f4", "_FillValue": None},
                            "time": {"dtype": "f4", "_FillValue": None}})

    summary = run_checker(case_dir)

    assert summary["total_spatial_errors"] >= 1
    assert "nodes, expected" in summary["log_text"]


def test_wrong_lat_range(case_dir):
    """Replace a 3D file with one where lat spans [0, 180] instead of [-90, 90]."""
    path = file_for(case_dir, "delta_g")
    path.unlink()

    bad_lat = np.linspace(0.0, 180.0, 257, dtype=np.float32)
    lon = np.linspace(0.0, 360.0, 513, endpoint=False, dtype=np.float32)
    time_years = np.array([START_YEAR, END_YEAR], dtype=np.float32)
    data = np.random.uniform(-10.0, 10.0, (2, 257, 513)).astype(np.float32)

    ds = _make_spatial_dataset("delta_g", data, time_years, bad_lat, lon)
    ds.to_netcdf(path, unlimited_dims=("time",),
                 encoding={"delta_g": {"dtype": "f4", "_FillValue": None},
                            "time": {"dtype": "f4", "_FillValue": None}})

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
        workdir=str(REPO_ROOT),
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
        workdir=str(REPO_ROOT),
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
    del_global_attr(first_file(case_dir), "contact_email")

    summary = run_checker(case_dir)

    assert summary["total_meta_attr_errors"] >= 1
    assert "contact_email" in summary["log_text"]


def test_wrong_reference_frame(case_dir):
    set_global_attr(first_file(case_dir), "reference_frame", "CE")

    summary = run_checker(case_dir)

    assert summary["total_meta_attr_errors"] >= 1
    assert "reference_frame" in summary["log_text"]


def test_missing_long_name(case_dir):
    path = file_for(case_dir, "maf")
    del_var_attr(path, "maf", "long_name")

    summary = run_checker(case_dir)

    assert summary["total_var_attr_errors"] >= 1
    assert "long_name" in summary["log_text"]


def test_wrong_dtype(case_dir, tmp_path):
    """Replace the delta_bed file with one where the variable is float64."""
    path = file_for(case_dir, "delta_bed")
    path.unlink()

    lat = np.linspace(-90.0, 90.0, 257, dtype=np.float32)
    lon = np.linspace(0.0, 360.0, 513, endpoint=False, dtype=np.float32)
    time_years = np.array([START_YEAR, END_YEAR], dtype=np.float32)
    data = np.random.uniform(-100.0, 100.0, (2, 257, 513)).astype(np.float64)

    ds = _make_spatial_dataset("delta_bed", data, time_years, lat, lon)
    ds.to_netcdf(path, unlimited_dims=("time",),
                 encoding={"delta_bed": {"dtype": "f8", "_FillValue": None},
                            "time": {"dtype": "f4", "_FillValue": None}})

    summary = run_checker(case_dir)

    assert summary["total_var_attr_errors"] >= 1
    assert "float32" in summary["log_text"]


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

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
            "time": ("time", time_years, {"units": generator.TIME_UNITS}),
            "lat": ("lat", lat, {"units": "degrees_north"}),
            "lon": ("lon", lon, {"units": "degrees_east"}),
        },
    )
    ds.attrs.update({
        "group": GROUP, "model": MODEL,
        "contact_name": "Synthetic Test", "contact_email": "test@example.org",
        "reference_frame": "CM",
    })
    return ds
