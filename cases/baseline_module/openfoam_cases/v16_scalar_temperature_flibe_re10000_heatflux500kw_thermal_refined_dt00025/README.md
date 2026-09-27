# v16 Scalar Temperature FLiBe Re10000 Heat Flux 500 kW Thermal Refined dt00025

High heat-flux passive-temperature sensitivity on the second-refinement mesh.

## Purpose

Repeat the thermal-refined FLiBe-like passive scalar case with five times the wall heat flux used in `v14`.

This checks whether the current OpenFOAM workflow remains stable and interpretable under a more aggressive blanket-relevant heat load.

## Physics

- clean setup copied from `v14` `0`, `constant`, and `system`
- frozen velocity/flux field from `v13`
- scalar field `T`, interpreted as temperature in K
- inlet T = 800 K
- FLiBe-like placeholder thermal diffusivity: 2.1478e-7 m^2/s
- heated-wall fixed gradient equivalent to about 500 kW/m^2 with k = 1.0 W/m/K
- deltaT = 0.00025 s
- target final time = 2 s

## Expected Use

Compare against `v14`:

- outlet average T
- heated-wall average T
- min/max T
- wall-surface average and max T - Tin
- line-profile shape and max T - Tin

## Current Result

Run completed to t = 2 s.

Final post-processing:

- outlet average T: about 801.022 K
- wall-surface average T: about 1018.888 K
- wall-surface max T: about 1116.608 K
- min T: about 799.552 K
- max T: about 1116.608 K

Generated outputs:

- `../../../../results/visualization/v16_patch_temperature_summary.csv`
- `../../../../results/heat_flux_sensitivity/flibe_thermal_refined_heat_flux_summary.csv`
- `../../../../results/heat_flux_sensitivity/flibe_thermal_refined_heat_flux_summary.png`
- `../../../../results/heat_flux_sensitivity/comparison_flibe_thermal_refined_100kw_500kw_wall_deltaT.png`
- `../../../../results/wall_profiles/v16_flibe_thermal_refined_500kw_heated_wall_surface_profile.csv`
- `../../../../results/wall_profiles/v16_flibe_thermal_refined_500kw_heated_wall_surface_profile.png`
- `../../../../results/line_profiles/v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025_vertical_chamber_midplane.csv`
- `../../../../results/line_profiles/v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025_vertical_chamber_midplane.png`
- `../../../../results/line_profiles/v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025_streamwise_near_heated_wall.csv`
- `../../../../results/line_profiles/v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025_streamwise_near_heated_wall.png`
- `VTK/v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025_8000.vtm`

Interpretation:

The 500 kW/m^2 response is close to linear relative to `v14`. Wall-surface average T - Tin increased from about 43.8 K to about 218.9 K, and outlet average T - Tin increased from about 0.20 K to about 1.02 K.
