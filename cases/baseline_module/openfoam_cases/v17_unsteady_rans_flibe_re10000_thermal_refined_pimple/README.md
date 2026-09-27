# v17 Unsteady RANS FLiBe Re10000 Thermal Refined PIMPLE

First transient-flow case on the current thermal-refined geometry.

## Purpose

This is the conservative first step toward the planned RANS-plus-unsteady comparison. It is not LES yet. It uses the converged `v13` steady RANS field as the initial condition and switches the flow solver from `simpleFoam` to `pimpleFoam`.

## Case Setup

- solver: `pimpleFoam`
- turbulence model: k-omega SST RANS
- mesh: copied from `v13_isothermal_rans_flibe_re10000_thermal_refined`
- initial fields: copied from the converged `v13` time `1000`
- time scheme: first-order Euler
- initial smoke-test duration: `0.02 s`
- adaptive time step: `maxCo 0.5`, `maxDeltaT 2.5e-4 s`
- post-processing dictionary: `system/fieldAverages`

## Smoke-Test Result

The short run completed successfully.

- latest saved time: `0.0201038 s`
- saved transient fields: `0.0101038`, `0.0151038`, `0.0201038`
- wall-clock runtime: about `23 s`
- final Courant number: mean about `0.074`, max about `0.474`
- final inlet areaAverage(p): `0.916451 m^2/s^2`
- final outlet areaAverage(p): `0 m^2/s^2`
- pressure-drop estimate using rho = 1940 kg/m^3: about `1778 Pa`
- final inlet flow rate: `-0.00092784 m^3/s`
- final outlet flow rate: `0.00092784 m^3/s`

This is close to the `v13` steady pressure-drop estimate of about `1774 Pa`, which is expected because the transient starts from the converged steady field and only runs for a short time.

Logs:

- `log.pimpleFoam`
- `log.postProcess_fieldAverages`
- `log.postProcess_cellCentres`
- `log.foamToVTK`

Generated outputs:

- `results/slices/v17_unsteady_smoke_midplane_z_velocity.png`
- `results/slices/v17_unsteady_smoke_midplane_z_velocity.csv`
- `VTK/v17_unsteady_rans_flibe_re10000_thermal_refined_pimple_82.vtm`

## Why This Exists

The current project already has a steady RANS flow baseline plus passive thermal scalar cases. Before attempting LES, this case checks whether the same geometry, boundary conditions, and turbulence fields can run as a stable transient calculation.

If this case starts cleanly, the next useful version is a longer unsteady RANS run long enough to compare recirculation and hot-region persistence against the steady baseline. LES should wait until this transient case has a known runtime and stable residual behavior.

## What To Compare

Compare against:

- `v13_isothermal_rans_flibe_re10000_thermal_refined` for steady velocity and pressure drop
- `v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025` for passive thermal response at 100 kW/m^2
- `v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025` for heat-flux sensitivity context

Useful transient checks:

- residual stability
- Courant number staying under the target
- outlet flow rate remaining close to the steady baseline
- whether recirculation structures remain steady or begin oscillating
- whether a later passive-scalar replay needs a time-dependent velocity field

## Claim Boundary

This case can support: "a transient RANS case was initialized from the steady baseline as the first step toward unsteady comparison."

This case cannot support: "LES was completed," "the flow is experimentally validated," or "the wall-temperature field is grid independent."
