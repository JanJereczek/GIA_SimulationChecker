# GIAMIP Variable Reference

Column definitions follow the GIAMIP data request spreadsheet. The **Precision** column reflects this checker's convention: all variables are stored as **float32** (single precision). The original GIAMIP specification requests float64 for most variables; if your group follows the original specification, please adjust accordingly.

---

## Required variables

| Description | Dimensions | Variable name | `standard_name` | `long_name` | Output interval | Units | Precision | `_FillValue` | Reference |
|-------------|-----------|---------------|-----------------|-------------|-----------------|-------|-----------|--------------|-----------|
| Bed / sea floor | time, lat, lon | `bed` | `bedrock_altitude` | Height of the solid Earth surface beneath ice and ocean water | 1000 years | m | float32 | none | Present-day reference ellipsoid / center of mass |
| Changes in geoid height | time, lat, lon | `delta_g` | — | Change in geoid height relative to the initial simulation time step | 1000 years | m | float32 | none | Geoid at initial timestep |
| Changes in relative sea level | time, lat, lon | `delta_rsl` | `change_in_mean_sea_level_wrt_solid_surface` | Relative sea-level change relative to the initial simulation timestep | 1000 years | m | float32 | none | Initial simulation timestep |
| Ocean area fraction | time, lat, lon | `ocean_area_fraction` | `sea_area_fraction` | Fraction of horizontal grid-cell area covered by ocean | Forcing interval | 1 | float32 | none | — |
| Total ice area fraction | time, lat, lon | `land_ice_area_fraction` | `land_ice_area_fraction` | Fraction of grid-cell area covered by grounded and floating land ice | Forcing interval | 1 | float32 | none | — |
| Spatial mean geoid change over ocean | time | `mean_delta_g` | — | Spatial mean of geoid height change over the ocean area | 1000 years | m | float32 | none | Initial timestep geoid |
| Grounded ice mass | time | `grd_ice_mass` | — | Spatial integration of grounded ice volume times ice density | Forcing interval | kg | float32 | none | — |
| Total ice mass | time | `total_ice_mass` | — | Total ice volume times ice density | Forcing interval | kg | float32 | none | — |
| Total ocean area (incl. grounded ice) | time | `ocean_area_grdice` | — | Total ocean area including marine regions covered by grounded ice | Forcing interval | m² | float32 | none | — |
| Total ocean area (excl. grounded ice) | time | `ocean_area` | — | Total ocean area excluding marine regions covered by grounded ice | Forcing interval | m² | float32 | none | — |
| Mass above flotation | time | `maf` | `land_ice_mass_not_displacing_sea_water` | Land ice mass above flotation | Forcing interval | kg | float32 | none | — |

---

## Optional variables

| Description | Dimensions | Variable name | `standard_name` | `long_name` | Output interval | Units | Precision | `_FillValue` | Reference |
|-------------|-----------|---------------|-----------------|-------------|-----------------|-------|-----------|--------------|-----------|
| Eastward horizontal solid Earth displacement | time, lat, lon | `delta_bed_east` | — | Eastward horizontal solid Earth displacement | 1000 years | m | float32 | none | Initial simulation timestep |
| Northward horizontal solid Earth displacement | time, lat, lon | `delta_bed_north` | — | Northward horizontal solid Earth displacement | 1000 years | m | float32 | none | Initial simulation timestep |
| Cosine spherical harmonic coefficients | degree, order | `Clm` | — | Cosine spherical harmonic coefficients of geoid height change | Once (snapshot) | 1 | float32 | none | — |
| Sine spherical harmonic coefficients | degree, order | `Slm` | — | Sine spherical harmonic coefficients of geoid height change | Once (snapshot) | 1 | float32 | none | — |

---

## Grid specification

| Property | Value |
|----------|-------|
| Grid type | Global Gaussian lat/lon |
| Latitude nodes | 257 |
| Longitude nodes | 513 |
| Latitude range | −90° to +90° (south to north) |
| Longitude range | 0° to <360° (west to east) |
| Time encoding | `days since 0001-01-01 00:00:00`, calendar `proleptic_gregorian` |
| Reference frame | CM (center of mass) |

---

## Notes

- **No missing values**: all data arrays must be free of NaN and fill values.
- **Mask variables** (`ocean_area_fraction`, `land_ice_area_fraction`): values must lie within [0, 1].
- **Spherical harmonics** (`Clm`, `Slm`): stored as an upper-triangular (degree, order) matrix of size 97 × 97 (degrees 0–96). No time or spatial dimensions.
- **Time step tolerance** for 1000-year variables: steps in the range [900, 1100] years are accepted.
