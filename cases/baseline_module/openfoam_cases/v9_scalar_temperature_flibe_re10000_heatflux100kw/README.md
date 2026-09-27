# v9 Scalar Temperature FLiBe Re10000 Heat Flux 100 kW

First FLiBe-like passive-temperature transport case.

## Purpose

Use the converged FLiBe-like Re = 10,000 velocity field from `v8` and solve passive scalar `T`.

This is still not buoyant flow, MHD, or conjugate heat transfer. It is the first molten-salt-like thermal transport debug case.

## Physics

- frozen flow field from `v8`
- scalar field `T`, interpreted as temperature in K
- FLiBe-like thermal diffusivity: `2.1478e-7 m^2/s`
- inlet scalar: 800 K
- heated wall: fixed gradient equivalent to about 100 kW/m^2 with thermal conductivity 1.0 W/m/K
- adiabatic walls: zero gradient
- outlet: inletOutlet

## Visualization

After running, open this in ParaView:

- `VTK/v9_scalar_temperature_flibe_re10000_heatflux100kw_1000.vtm`

## Current Result

Run completed to t = 2 s.

Final post-processing:

- outlet average T: about 800.123 K
- heated-wall average T: about 1084.43 K
- min T: 799.669 K
- max T: 1119.34 K

Generated line-profile outputs:

- `../../../../results/line_profiles/v9_scalar_temperature_flibe_re10000_heatflux100kw_vertical_chamber_midplane.csv`
- `../../../../results/line_profiles/v9_scalar_temperature_flibe_re10000_heatflux100kw_vertical_chamber_midplane.png`
- `../../../../results/line_profiles/v9_scalar_temperature_flibe_re10000_heatflux100kw_streamwise_near_heated_wall.csv`
- `../../../../results/line_profiles/v9_scalar_temperature_flibe_re10000_heatflux100kw_streamwise_near_heated_wall.png`

Important interpretation note:

The line profiles are nearest-cell interior samples. They are useful for early comparison plots, but they do not replace wall-surface temperature extraction or a near-wall mesh refinement check.
