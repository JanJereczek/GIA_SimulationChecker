# GIAMIP GIA Simulation Compliance Checker

Checks GIAMIP NetCDF simulation datasets for compliance with the GIAMIP data request conventions. Each experiment is first checked for **missing mandatory variables**, then the following categories are validated for every file:

1. **Naming** — variable name, experiment ID, group name, and model name follow the filename convention.
2. **Numerical** — no NaN/missing values; all values lie within the variable's plausible range (`bounds` in the variable metadata).
3. **Spatial** *(lat/lon/time variables only)* — global Gaussian grid of 257 × 513 nodes; latitude spans [−90, 90]; longitude spans [0, 360); and the lat/lon grid matches the forcing grid (within a relative tolerance).
   - **Dimension ordering** — variable dimensions are in the expected order: `(time, lat, lon)` for fields, `(time,)` for scalars, `(degree, order)` for `Clm`/`Slm`.
4. **Time** — time axis is monotonically increasing; time values are calendar years and are compared year-by-year against the forcing file's `year` variable.
5. **Attributes** — reported in three buckets: **metadata** (global attributes `group`, `model`, `contact_name`, `contact_email`, `reference_frame`; `reference_frame` must be `CM`), **coordinate** (time units `year`; lat/lon units present), and **variable** (units match the data request, `long_name` present, `standard_name` matches when specified, data variable is float32).

Compliance criteria are defined in `giamip_compliance_checker.py`: variable metadata in `GIAMIP_VARIABLES` and the valid experiments (`Exp01`–`Exp12`) in `GIAMIP_EXPERIMENTS`. The forcing file supplies the reference time axis and grid.

---

## Setup

To create the environment (only the first time):
```bash
conda env create -f isschecker_env.yml
```

Once the environment is created, you need to activate it (each time you use the code):
```bash
conda activate isschecker
```

Dependencies: Python 3.14, `numpy` 2.4, `xarray` 2026.4, `netCDF4` 1.7, `tqdm` 4.67.

---

## Running the checker

All scripts must be run from the repository root. There are three entry points, at increasing scope. Argument names indicate whether a **directory** (`-dir`) or a single **file** (`-filepath`) is expected.

### 1. Single experiment

Checks the `.nc` files in **one experiment directory** against **one forcing file**. Writes `compliance_checker_log.txt` into the source directory.

```bash
python giamip_compliance_checker.py \
    --source-dir ./output/CUBoulder-SemiAnalytic/Exp05 \
    --forcing-filepath ./input/iceHistory-PaleoMIST_1a.nc
```

| Option | Default | Description |
|--------|---------|-------------|
| `--source-dir` | *(required)* | Directory of `.nc` files for a single experiment |
| `--forcing-filepath` | *(required)* | Forcing NetCDF **file** used to validate output time axes and grid |

### 2. All experiments of one model

Loops over every experiment subdirectory of a **model directory** (`Exp01`, `Exp05`, …) and automatically selects each experiment's forcing file from `--forcing-dir`, according to the mapping in [experiments.md](experiments.md). Missing forcing files are skipped with a warning (the experiment is still checked, minus the time/grid comparisons).

```bash
python giamip_model_compliance_checker.py \
    --source-dir ./output/CUBoulder-SemiAnalytic \
    --forcing-dir ./input/
```

| Option | Default | Description |
|--------|---------|-------------|
| `--source-dir` | *(required)* | Model directory containing experiment subdirectories |
| `--forcing-dir` | *(required)* | Directory containing the forcing NetCDF files |

### 3. All models at once

Loops over every model directory under the output directory, each looping over its experiments (as above). Use this to check all submitted experiments in one go.

```bash
python giamip_full_compliance_checker.py \
    --source-dir ./output \
    --forcing-dir ./input
```

| Option | Default | Description |
|--------|---------|-------------|
| `--source-dir` | *(required)* | Top-level directory containing model subdirectories |
| `--forcing-dir` | *(required)* | Directory containing the forcing NetCDF files |

---

## File naming convention

Each output file must be named:

```
{variable_name}_{experiment_id}_{group}_{model}.nc
```

For example: `delta_bed_Exp01_AWI_MyModel.nc`

Valid experiment IDs (`Exp01`–`Exp12`) are defined by `GIAMIP_EXPERIMENTS` in `giamip_compliance_checker.py`.

---

## Variables

Full variable metadata (long names, standard names, units, precision, reference surface) is in [variables.md](variables.md).

| Variable | Dimensions | Mandatory | Output interval |
|----------|-----------|-----------|-----------------|
| `delta_bed` | lat, lon, time | yes | forcing |
| `delta_g` | lat, lon, time | yes | forcing |
| `ocean_area_fraction` | lat, lon, time | yes | forcing |
| `land_ice_area_fraction` | lat, lon, time | yes | forcing |
| `mean_delta_g` | time | yes | forcing |
| `grd_ice_mass` | time | yes | forcing |
| `total_ice_mass` | time | yes | forcing |
| `ocean_area_grdice` | time | yes | forcing |
| `ocean_area` | time | yes | forcing |
| `maf` | time | yes | forcing |
| `delta_bed_east` | lat, lon, time | no | forcing |
| `delta_bed_north` | lat, lon, time | no | forcing |
| `Clm` | degree, order | no | snapshot |
| `Slm` | degree, order | no | snapshot |

---

## Generating synthetic test files

`generate/generate_giamip_test_files.py` creates GIAMIP-compliant NetCDF test files with synthetic data.

```bash
conda activate isschecker

# Generate all mandatory variables for Exp01
python generate/generate_giamip_test_files.py --experiment-id Exp01 --group MYGROUP --model MyModel

# Generate a single variable
python generate/generate_giamip_test_files.py --variable delta_bed --experiment-id Exp01 --group MYGROUP --model MyModel

# Custom time range and step count
python generate/generate_giamip_test_files.py --start-year 1 --end-year 2001 --n-steps 3
```

| Option | Default | Description |
|--------|---------|-------------|
| `--output-dir` | `./output/GIAMIP/Exp01/CORE` | Directory for output files |
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
