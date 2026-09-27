# v15 Scalar Temperature FLiBe Re10000 Heat Flux 100 kW Thermal Refined Long Run dt00025

Longer passive-temperature continuation of the second-refinement FLiBe-like 100 kW/m^2 case.

## Purpose

Extend `v14` past t = 2 s to check whether outlet temperature and heated-wall temperature are still evolving over the short passive-scalar window.

This is not a new mesh. It is a time-extension check using the same mesh and frozen flow as `v14`.

## Physics

- copied from completed `v14` result
- frozen velocity/flux field from `v13`
- scalar field `T`, interpreted as temperature in K
- inlet T = 800 K
- FLiBe-like placeholder thermal diffusivity: 2.1478e-7 m^2/s
- heated-wall fixed gradient equivalent to about 100 kW/m^2 with k = 1.0 W/m/K
- deltaT = 0.00025 s
- starts from latest available time, initially t = 2 s
- target final time = 5 s

## Expected Use

Use this case to decide whether a later 10-20 s run is worth the compute time.

Key quantities to compare against `v14` at t = 2 s:

- outlet average T
- heated-wall average T
- min/max T
- wall-surface average and max T - Tin
- shape of the near-wall streamwise temperature profile

## Current Result

The solver marched to t = 5.0 s. Because the write interval lands on 4.8 s, the latest saved field used for post-processing is t = 4.8 s.

Final saved-field post-processing:

- outlet average T: 800.213 K
- heated-wall average T: 843.848 K
- min T: 799.910 K
- max T: 863.595 K
- wall surface average T - Tin: 43.85 K
- wall surface max T - Tin: 63.59 K

VTK:

- `VTK/v15_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_longrun_dt00025_19200.vtm`

Generated outputs:

- `../../../../results/visualization/v15_patch_temperature_summary.csv`
- `../../../../results/wall_profiles/v15_flibe_thermal_refined_100kw_longrun_t4p8_heated_wall_surface_profile.csv`
- `../../../../results/wall_profiles/v15_flibe_thermal_refined_100kw_longrun_t4p8_heated_wall_surface_profile.png`
- `../../../../results/line_profiles/v15_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_longrun_dt00025_vertical_chamber_midplane.csv`
- `../../../../results/line_profiles/v15_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_longrun_dt00025_vertical_chamber_midplane.png`
- `../../../../results/line_profiles/v15_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_longrun_dt00025_streamwise_near_heated_wall.csv`
- `../../../../results/line_profiles/v15_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_longrun_dt00025_streamwise_near_heated_wall.png`

Interpretation:

The field changed very little relative to `v14` at t = 2 s. This makes the 500 kW/m^2 sensitivity a better next run than immediately spending more compute on a 10-20 s passive continuation.
