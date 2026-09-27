# v11 Scalar Temperature FLiBe Re10000 Heat Flux 100 kW Refined

Refined-mesh passive-temperature attempt using the same time step as coarse `v9`.

## Purpose

Test whether the coarse thermal setup could run unchanged on the refined mesh.

## Physics

- frozen refined velocity/flux field from `v10`
- scalar field `T`, interpreted as temperature in K
- inlet T = 800 K
- heated-wall fixed gradient equivalent to about 100 kW/m^2
- thermal diffusivity: 2.1478e-7 m^2/s
- deltaT = 0.002 s

## Outcome

This case did not complete.

Observed behavior:

- initial max Courant number: about 2.16
- solver reached about t = 0.774 s
- run ended with a floating-point exception in the scalar solve

Interpretation:

This is not a valid thermal result. It is useful evidence that the refined mesh needs a smaller scalar-transport time step. The stable repeat is `v12_scalar_temperature_flibe_re10000_heatflux100kw_refined_dt0005`.
