# v14 Scalar Temperature FLiBe Re10000 Heat Flux 100 kW Thermal Refined dt00025

Second-refinement passive-temperature repeat of the FLiBe-like 100 kW/m^2 case.

## Purpose

Repeat the refined thermal case after increasing heated-wall resolution again. This checks whether wall temperature is moving toward a stable value after coarse `v9` and refined `v12`.

## Physics

- frozen velocity/flux field from `v13`
- scalar field `T`, interpreted as temperature in K
- inlet T = 800 K
- FLiBe-like placeholder thermal diffusivity: 2.1478e-7 m^2/s
- heated-wall fixed gradient equivalent to about 100 kW/m^2 with k = 1.0 W/m/K
- deltaT = 0.00025 s
- final time = 2 s

## Current Result

Run completed to t = 2 s.

Final post-processing:

- outlet average T: about 800.200 K
- heated-wall average T: about 843.798 K
- min T: 799.910 K
- max T: 863.519 K

VTK:

- `VTK/v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025_8000.vtm`

Generated line-profile outputs:

- `../../../../results/line_profiles/v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025_vertical_chamber_midplane.csv`
- `../../../../results/line_profiles/v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025_vertical_chamber_midplane.png`
- `../../../../results/line_profiles/v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025_streamwise_near_heated_wall.csv`
- `../../../../results/line_profiles/v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025_streamwise_near_heated_wall.png`

Wall-profile outputs:

- `../../../../results/wall_profiles/v14_flibe_thermal_refined_100kw_heated_wall_surface_profile.csv`
- `../../../../results/wall_profiles/v14_flibe_thermal_refined_100kw_heated_wall_surface_profile.png`
- `../../../../results/wall_profiles/comparison_flibe_heated_wall_refinement_deltaT.png`

Interpretation:

The heated-wall average temperature dropped again relative to `v12`, from about 869 K to about 844 K. Wall temperature is still mesh-sensitive, but it is moving toward a much cooler value than coarse `v9`.
