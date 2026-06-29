# GIAMIP Variable Reference

Column definitions follow the GIAMIP data request spreadsheet (Version 2, 2026-06-08). The **Precision** column reflects this checker's convention: all variables are stored as **float32** (single precision). The original GIAMIP specification requests float64 for scalar time series and spherical harmonic coefficients; if your group follows the original specification, please adjust accordingly.

---

## Required variables

| Description | Dimensions | Variable name | `standard_name` | `long_name` | Output interval | Units | Plausible range | Range rationale | Precision | `_FillValue` | Reference |
|-------------|-----------|---------------|-----------------|-------------|-----------------|-------|-----------------|-----------------|-----------|--------------|-----------|
| Changes in bedrock elevation | time, lat, lon | `delta_bed` | — | Change in the bedrock elevation relative to the initial simulation time step | Forcing interval | m | [−3000, 3000] | Vertical GIA motion; LGM extremes ~1–2 km, so ±2 km is comfortably loose. | float32 | none | Bedrock elevation at initial simulation timestep |
| Changes in geoid height | time, lat, lon | `delta_g` | — | Change in the geoid height relative to the initial simulation time step | Forcing interval | m | [−500, 500] | Geoid change is smaller than bedrock motion; tens of m typical, up to hundreds locally. | float32 | none | Geoid at initial simulation timestep |
| Ocean area fraction | time, lat, lon | `ocean_area_fraction` | `sea_area_fraction` | Fraction of horizontal grid-cell area covered by ocean | Forcing interval | 1 | [-1e-6, 1 + 1e-6] | Fraction; already enforced as a mask variable. | float32 | none | — |
| Total ice area fraction | time, lat, lon | `land_ice_area_fraction` | `land_ice_area_fraction` | Fraction of horizontal grid-cell area covered by grounded and floating land ice | Forcing interval | 1 | [-1e-6, 1 + 1e-6] | Fraction; already enforced as a mask variable. | float32 | none | — |
| Spatial mean geoid change over ocean | time | `mean_delta_g` | — | Spatial mean of geoid height change (delta_g) over the ocean area | Forcing interval | m | [−300, 300] | Ocean-mean (barystatic) geoid change is O(100 m) over a glacial cycle. | float32 | none | Geoid at initial simulation timestep |
| Grounded ice mass | time | `grd_ice_mass` | — | Spatial integration of grounded ice volume times ice density | Forcing interval | kg | [1e18, 1e20] | Mass ≥ 0; present grounded ice ≈ 2.7e19 kg, LGM ~2–3×. | float32 | none | — |
| Total ice mass | time | `total_ice_mass` | — | Spatial integration, total (grounded and floating) ice volume times ice density | Forcing interval | kg | [1e18, 1e20] | Slightly above grounded ice (adds floating); same generous ceiling. | float32 | none | — |
| Total ocean area (incl. grounded ice) | time | `ocean_area_grdice` | — | Total ocean area including marine regions covered by grounded ice | Forcing interval | m² | [1e14, 5e14] | ≤ Earth's surface area (5.1e14 m²); ocean ≈ 3.6e14 m². | float32 | none | — |
| Total ocean area (excl. grounded ice) | time | `ocean_area` | — | Total ocean area excluding marine regions covered by grounded ice | Forcing interval | m² | [1e14, 5e14] | Same physical ceiling as `ocean_area_grdice`. | float32 | none | — |
| Mass above flotation | time | `maf` | `land_ice_mass_not_displacing_sea_water` | Land ice mass above flotation that would contribute to global mean sea-level change if converted to water and added to the ocean | Forcing interval | kg | [1e19, 1e20] | Mass above flotation ≥ 0; order of total ice mass. | float32 | none | — |

---

## Optional variables (strongly recommended)

| Description | Dimensions | Variable name | `standard_name` | `long_name` | Output interval | Units | Plausible range | Range rationale | Precision | `_FillValue` | Reference |
|-------------|-----------|---------------|-----------------|-------------|-----------------|-------|-----------------|-----------------|-----------|--------------|-----------|
| Eastward horizontal solid Earth displacement | time, lat, lon | `delta_bed_east` | — | Eastward horizontal solid Earth displacement relative to the initial simulation timestep | Forcing interval | m | [−200, 200] | Horizontal displacement ≪ vertical (realistically ≤ tens of m). | float32 | none | Geoid at initial simulation timestep |
| Northward horizontal solid Earth displacement | time, lat, lon | `delta_bed_north` | — | Northward horizontal solid Earth displacement relative to the initial simulation timestep | Forcing interval | m | [−200, 200] | Same as the eastward component. | float32 | none | Geoid at initial simulation timestep |
| Cosine Stokes coefficients of geoid change | degree, order | `Clm` | — | Cosine spherical harmonic coefficients (C_lm) of geoid height change (delta_g) between the first and final simulation timesteps | Once (snapshot) | 1 | — | No bound: observed magnitudes are O(100), so the normalization/units are not pinned down enough to set a meaningful range. | float32 | none | Initial simulation timestep |
| Sine Stokes coefficients of geoid change | degree, order | `Slm` | — | Sine spherical harmonic coefficients (S_lm) of geoid height change (delta_g) between the first and final simulation timesteps | Once (snapshot) | 1 | — | Same as `Clm`. | float32 | none | Initial simulation timestep |

---

## Grid specification

| Property | Value |
|----------|-------|
| Grid type | Global Gaussian lat/lon |
| Latitude nodes | 257 |
| Longitude nodes | 513 |
| Latitude range | −90° to +90° (south to north) |
| Longitude range | 0° to <360° (west to east) |
| Time encoding | years before year 2000 |
| Reference frame | CM (center of mass) |

---

## Notes

- **Plausible range**: loose sanity-check bounds intended to catch gross errors (wrong units, sign bugs, garbage fill values), not tight physical limits. Signed fields assume symmetric ± bounds; mass/area fields are non-negative. These are documentation only and are not yet enforced by the checker.
- **No missing values**: all data arrays must be free of NaN and fill values.
- **Mask variables** (`ocean_area_fraction`, `land_ice_area_fraction`): values must lie within [0, 1]; conservative remapping recommended for interpolation.
- **`total_ice_mass`**: should equal floating ice mass plus grounded ice mass.
- **`maf`**: barystatic sea-level change is derived from MAF divided by seawater density and ocean area.
- **Spherical harmonics** (`Clm`, `Slm`): stored as an upper-triangular (degree, order) matrix of size 97 × 97 (degrees 0–96), using the fully_normalized_4pi convention.
- **`delta_bed_east` / `delta_bed_north`**: groups with horizontal deformation capabilities are strongly encouraged to provide these.
- **Ocean depth change**: can be computed as `delta_bed − delta_g`.
