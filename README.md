# GIAMIP GIA Simulation Compliance Checker

Checks GIAMIP NetCDF simulation datasets for compliance with the GIAMIP data request conventions. The following categories are validated for every file:

1. **Naming** — variable name, experiment ID, group name, and model name; mandatory variables must all be present.
2. **Numerical** — units match the data request; no NaN values are present; mask variables (`ocean_area_fraction`, `land_ice_area_fraction`) have all values within [0, 1].
3. **Spatial** *(lat/lon/time variables only)* — global Gaussian grid of 257 × 513 nodes; latitude spans [−90, 90]; longitude spans [0, 360).
4. **Time** — time axis is monotonically increasing; start and end years fall within the allowed range for the experiment; variables with `1000yr` output interval use approximately 1000-year time steps.
5. **Attributes** — required global attributes (`group`, `model`, `contact_name`, `contact_email`, `reference_frame`) are present; `reference_frame` must be `CM`; time units start with `days since`; all data variables are float32.

Compliance criteria are defined in `giamip_compliance_checker.py` (variable metadata) and `experiments_giamip.csv` (valid experiment year ranges).

---

## Setup

```bash
conda env create -f isschecker_env.yml
conda activate isschecker
```

Dependencies: Python 3.14, `numpy` 2.4, `xarray` 2026.4, `cftime` 1.6, `netCDF4` 1.7, `tqdm` 4.67.

---

## Running the checker

The script must be run from the repository root. It writes `giamip_compliance_checker_log.txt` into the `--source-path` directory.

```bash
python giamip_compliance_checker.py --source-path ./Models/GIAMIP/Exp01/CORE
```

| Option | Default | Description |
|--------|---------|-------------|
| `--source-path` | `./Models/GIAMIP/Exp01/CORE` | Directory containing `.nc` files to check |

---

## File naming convention

Each output file must be named:

```
{variable_name}_{experiment_id}_{group}_{model}.nc
```

For example: `bed_Exp01_AWI_MyModel.nc`

Valid experiment IDs are defined in `experiments_giamip.csv`.

---

## Variables

Full variable metadata (long names, standard names, units, precision, reference surface) is in [variables.md](variables.md).

| Variable | Dimensions | Mandatory | Output interval |
|----------|-----------|-----------|-----------------|
| `bed` | lat, lon, time | yes | 1000 yr |
| `maf` | lat, lon, time | yes | 1000 yr |
| `rsl` | lat, lon, time | yes | 1000 yr |
| `delta_g` | lat, lon, time | yes | 1000 yr |
| `delta_rsl` | lat, lon, time | yes | 1000 yr |
| `ocean_area_fraction` | lat, lon, time | yes | 1000 yr |
| `land_ice_area_fraction` | lat, lon, time | yes | 1000 yr |
| `grd_ice_mass` | time | yes | 1000 yr |
| `mean_delta_g` | time | yes | 1000 yr |
| `Clm` | degree, order | yes | snapshot |
| `Slm` | degree, order | yes | snapshot |
| `oaf_10yr` | lat, lon, time | no | 10 yr |
| `liaf_10yr` | lat, lon, time | no | 10 yr |
| `grd_ice_mass_10yr` | time | no | 10 yr |
| `mean_delta_g_10yr` | time | no | 10 yr |

---

## Generating synthetic test files

`generate/generate_giamip_test_files.py` creates GIAMIP-compliant NetCDF test files with synthetic data.

```bash
conda activate isschecker

# Generate all mandatory variables for Exp01
python generate/generate_giamip_test_files.py --experiment-id Exp01 --group AWI --model MyModel

# Generate a single variable
python generate/generate_giamip_test_files.py --variable bed --experiment-id Exp01 --group AWI --model MyModel

# Custom time range and step count
python generate/generate_giamip_test_files.py --start-year 1 --end-year 2001 --n-steps 3
```

| Option | Default | Description |
|--------|---------|-------------|
| `--output-dir` | `./Models/GIAMIP/Exp01/CORE` | Directory for output files |
| `--experiment-id` | `Exp01` | Experiment identifier |
| `--group` | `TESTGROUP` | Group name |
| `--model` | `TESTMODEL` | Model name |
| `--start-year` | `1` | First year of the time axis |
| `--end-year` | `1001` | Last year of the time axis |
| `--n-steps` | `2` | Number of time steps |
| `--variable` | *(all mandatory)* | Generate a single variable only |

---

## Running tests

The regression suite uses `pytest` and creates temporary synthetic datasets, then mutates them to verify expected checker failures for naming, missing variables, numerical, spatial, time-axis, and attribute errors.

```bash
pytest -v tests/test_giamip_compliance_checker.py
```

If you want to retain the files generated during testing you can use:

```bash
pytest -v tests/test_giamip_compliance_checker.py --basetemp=/tmp/pytest_tmp
```

The files will then be left in `/tmp/pytest_tmp`. Otherwise, they are cleaned up once tests pass.
