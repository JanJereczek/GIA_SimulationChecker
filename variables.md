# GIAMIP Variable Reference

Column definitions follow the GIAMIP data request spreadsheet (Version 3, 2026-07-02). The **Precision** column matches the data request: **float32** for the gridded (time, lat, lon) fields, **float64** for the scalar time series and the spherical harmonic coefficients. The checker enforces this per variable.

---

## Required variables

| Description | Dimensions | Variable name | `standard_name` | `long_name` | Output interval | Units | Plausible range | Range rationale | Precision | `_FillValue` | Reference |
|-------------|-----------|---------------|-----------------|-------------|-----------------|-------|-----------------|-----------------|-----------|--------------|-----------|
| Changes in bedrock elevation | time, lat, lon | `delta_bed` | — | Change in the bedrock elevation relative to the initial simulation time step | Forcing interval | m | [−3000, 3000] | Vertical GIA motion; LGM extremes ~1–2 km, so ±2 km is comfortably loose. | float32 | none | Bedrock elevation at initial simulation timestep |
| Changes in geoid height | time, lat, lon | `delta_g` | — | Change in the geoid height relative to the initial simulation time step | Forcing interval | m | [−500, 500] | Geoid change is smaller than bedrock motion; tens of m typical, up to hundreds locally. | float32 | none | Geoid at initial simulation timestep |
| Ocean area fraction | time, lat, lon | `ocean_area_fraction` | `sea_area_fraction` | Fraction of horizontal grid-cell area covered by ocean | Forcing interval | 1 | [-1e-6, 1 + 1e-6] | Fraction; already enforced as a mask variable. | float32 | none | — |
| Total ice area fraction | time, lat, lon | `land_ice_area_fraction` | `land_ice_area_fraction` | Fraction of horizontal grid-cell area covered by grounded and floating land ice | Forcing interval | 1 | [-1e-6, 1 + 1e-6] | Fraction; already enforced as a mask variable. | float32 | none | — |
| Spatial mean geoid change over ocean | time | `mean_delta_g` | — | Spatial mean of geoid height change (delta_g) over the ocean area | Forcing interval | m | [−300, 300] | Ocean-mean (barystatic) geoid change is O(100 m) over a glacial cycle. | float64 | none | Geoid at initial simulation timestep |
| Grounded ice mass | time | `grd_ice_mass` | — | Spatial integration of grounded ice volume times ice density | Forcing interval | kg | [1e18, 1e20] | Mass ≥ 0; present grounded ice ≈ 2.7e19 kg, LGM ~2–3×. | float64 | none | — |
| Total ice mass | time | `total_ice_mass` | — | Spatial integration, total (grounded and floating) ice volume times ice density | Forcing interval | kg | [1e18, 1e20] | Slightly above grounded ice (adds floating); same generous ceiling. | float64 | none | — |
| Total ocean area (incl. grounded ice) | time | `ocean_area_grdice` | — | Total ocean area including marine regions covered by grounded ice | Forcing interval | m² | [1e14, 5e14] | ≤ Earth's surface area (5.1e14 m²); ocean ≈ 3.6e14 m². | float64 | none | — |
| Total ocean area (excl. grounded ice) | time | `ocean_area` | — | Total ocean area excluding marine regions covered by grounded ice | Forcing interval | m² | [1e14, 5e14] | Same physical ceiling as `ocean_area_grdice`. | float64 | none | — |
| Mass above flotation | time | `maf` | `land_ice_mass_not_displacing_sea_water` | Land ice mass above flotation that would contribute to global mean sea-level change if converted to water and added to the ocean | Forcing interval | kg | [1e19, 1e20] | Mass above flotation ≥ 0; order of total ice mass. | float64 | none | — |

---

## Optional variables (strongly recommended)

| Description | Dimensions | Variable name | `standard_name` | `long_name` | Output interval | Units | Plausible range | Range rationale | Precision | `_FillValue` | Reference |
|-------------|-----------|---------------|-----------------|-------------|-----------------|-------|-----------------|-----------------|-----------|--------------|-----------|
| Eastward horizontal solid Earth displacement | time, lat, lon | `delta_bed_east` | — | Eastward horizontal solid Earth displacement relative to the initial simulation timestep | Forcing interval | m | [−200, 200] | Horizontal displacement ≪ vertical (realistically ≤ tens of m). | float32 | none | Geoid at initial simulation timestep |
| Northward horizontal solid Earth displacement | time, lat, lon | `delta_bed_north` | — | Northward horizontal solid Earth displacement relative to the initial simulation timestep | Forcing interval | m | [−200, 200] | Same as the eastward component. | float32 | none | Geoid at initial simulation timestep |
| Cosine Stokes coefficients of geoid change | degree, order | `Clm` | — | Cosine spherical harmonic coefficients (C_lm) of geoid height change (delta_g) between the first and final simulation timesteps | Once, at the last timestep | 1 | — | No bound: observed magnitudes are O(100), so the normalization/units are not pinned down enough to set a meaningful range. | float64 | none | Initial simulation timestep |
| Sine Stokes coefficients of geoid change | degree, order | `Slm` | — | Sine spherical harmonic coefficients (S_lm) of geoid height change (delta_g) between the first and final simulation timesteps | Once, at the last timestep | 1 | — | Same as `Clm`. | float64 | none | Initial simulation timestep |

---

## Coordinate variables

Version 3 of the data request adds an explicit list of coordinate variables, to be included in every output file. All of them carry `units` and `long_name`; those with an applicable CF name also carry `standard_name`.

| Description | Dimension | Variable name | `standard_name` | `long_name` | Units | Data type | Additional attributes |
|-------------|-----------|---------------|-----------------|-------------|-------|-----------|-----------------------|
| Time | time | `time` | `time` | Time | days since 2000-01-01 00:00:00 | float64 | `calendar = "365_day"` |
| Calendar year | time | `year` | — | Calendar year corresponding to each output time (positive = CE, negative = BCE) | year | int32 | — |
| Latitude | lat | `lat` | `latitude` | Latitude | degrees_north | float64 | — |
| Longitude | lon | `lon` | `longitude` | Longitude | degrees_east | float64 | — |
| Spherical harmonic degree | degree | `degree` | — | Spherical harmonic degree | 1 | float64 | — |
| Spherical harmonic order | order | `order` | — | Spherical harmonic order | 1 | float64 | — |

All of the above are checked, except that the presence of `degree`/`order` is not yet required (their data type is checked when they are there). `time` is the days axis and `year` carries the year of each step, offset from 2000 CE: `year == time/365`, so 0 -> 0, 365 -> 1, -365 -> -1. The checker verifies that correspondence, and uses `year` — not `time` — as the reference time vector when comparing against the forcing file.

CF treats `days since 2000-01-01` and `days since 2000-01-01 00:00:00` as identical, so both spellings are accepted.

---

## Global attributes

Every output file must carry the global attributes `group`, `model`, `contributors_name` and `contributors_email`. All four are checked.

---

## Grid specification

| Property | Value |
|----------|-------|
| Grid type | Global Gaussian lat/lon |
| Latitude nodes | 257 |
| Longitude nodes | 513 |
| Latitude range | −90° to +90° (south to north) |
| Longitude range | 0° to <360° (west to east) |
| Time encoding | `time` = days since 2000-01-01 (`365_day`); `year` = years offset from 2000 CE |

---

## Notes

- **Plausible range**: loose sanity-check bounds intended to catch gross errors (wrong units, sign bugs, garbage fill values), not tight physical limits. Signed fields assume symmetric ± bounds; mass/area fields are non-negative. These bounds are specific to this checker (they have no counterpart in the data request) and **are** enforced: values outside the range are reported as numerical errors.
- **`long_name` / `standard_name`**: the `long_name` of every data variable and of `time`/`lat`/`lon` must match the wording in the tables above (compared ignoring case and runs of whitespace, since the data request itself carries stray spaces). `standard_name` must match exactly, and is only checked for the variables that have one in the request — a `standard_name` declared on a variable whose entry is blank is not reported.
- **No missing values**: all data arrays must be free of NaN, infinities and fill values. The checker also rejects a declared `_FillValue` or `missing_value` attribute, since the data request specifies `_FillValue` "none" for every variable.
- **Mask variables** (`ocean_area_fraction`, `land_ice_area_fraction`): values must lie within [0, 1]; conservative remapping recommended for interpolation.
- **`total_ice_mass`**: should equal floating ice mass plus grounded ice mass.
- **`maf`**: barystatic sea-level change is derived from MAF divided by seawater density and ocean area.
- **Spherical harmonics** (`Clm`, `Slm`): stored as an upper-triangular (degree, order) matrix of size 97 × 97 (degrees 0–96), using the fully_normalized_4pi convention.
- **Dimension order in column-major languages**: the orders given above are the order the file stores, i.e. what `ncdump`, C and Python see. The netCDF Fortran and MATLAB interfaces reverse the dimension list, so declaring `(time, lat, lon)` there writes `(lon, lat, time)` to the file — declare them in reverse from those languages, and verify with `ncdump -h`. This trap is easy to miss for `Clm`/`Slm`, where both axes are 97 long and the sizes give nothing away.
- **`delta_bed_east` / `delta_bed_north`**: groups with horizontal deformation capabilities are strongly encouraged to provide these.
- **Ocean depth change**: can be computed as `delta_bed − delta_g`.
