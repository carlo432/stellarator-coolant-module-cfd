# v12 Scalar Temperature FLiBe Re10000 Heat Flux 100 kW Refined dt0005

Stable refined-mesh repeat of the FLiBe-like 100 kW/m^2 passive-temperature case.

## Purpose

Repeat coarse `v9` on the refined `v10` mesh with a smaller time step after `v11` showed scalar instability at `deltaT = 0.002 s`.

## Physics

- frozen refined velocity/flux field from `v10`
- scalar field `T`, interpreted as temperature in K
- inlet T = 800 K
- FLiBe-like placeholder thermal diffusivity: 2.1478e-7 m^2/s
- heated-wall fixed gradient equivalent to about 100 kW/m^2 with k = 1.0 W/m/K
- deltaT = 0.0005 s
- final time = 2 s

## Current Result

Run completed to t = 2 s.

Final post-processing:

- outlet average T: about 800.195 K
- heated-wall average T: about 869.171 K
- min T: 799.929 K
- max T: 893.437 K

VTK:

- `VTK/v12_scalar_temperature_flibe_re10000_heatflux100kw_refined_dt0005_4000.vtm`

Generated line-profile outputs:

- `../../../../results/line_profiles/v12_scalar_temperature_flibe_re10000_heatflux100kw_refined_dt0005_vertical_chamber_midplane.csv`
- `../../../../results/line_profiles/v12_scalar_temperature_flibe_re10000_heatflux100kw_refined_dt0005_vertical_chamber_midplane.png`
- `../../../../results/line_profiles/v12_scalar_temperature_flibe_re10000_heatflux100kw_refined_dt0005_streamwise_near_heated_wall.csv`
- `../../../../results/line_profiles/v12_scalar_temperature_flibe_re10000_heatflux100kw_refined_dt0005_streamwise_near_heated_wall.png`

Refinement comparison plots:

- `../../../../results/line_profiles/comparison_flibe_refinement_vertical_deltaT.png`
- `../../../../results/line_profiles/comparison_flibe_refinement_streamwise_deltaT.png`

Interpretation:

The refined thermal case is stable, but the thermal field is mesh-sensitive. Compared with coarse `v9`, the refined case has much lower wall average and maximum temperature at the same imposed heat flux and final time. Treat coarse wall-temperature values as provisional.
