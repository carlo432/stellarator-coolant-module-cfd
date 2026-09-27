# Baseline Case Matrix

Use this as the first simulation plan.

## Geometry

Coarse cases use:

- `cases/baseline_module/spec.md`
- `cases/baseline_module/geometry/baseline_module_v1.geo`

Refined cases use:

- `cases/baseline_module/spec.md`
- `cases/baseline_module/geometry/baseline_module_v2_refined.geo`
- `cases/baseline_module/geometry/baseline_module_v3_thermal_refined.geo`

## Target Cases

| Case | Fluid | Re | Heat Flux | Purpose |
|---|---:|---:|---:|---|
| A0 | water debug | 1,000 | 0 | mesh and boundary-condition debug |
| A1 | water debug | 10,000 | 0 | isothermal RANS baseline, done as `v5_isothermal_rans_water_re10000_clean` |
| A2 | water debug | 10,000 | 10 kW/m^2 | passive scalar thermal debug, done as `v6_scalar_temperature_water_re10000_heatflux` |
| A3 | water debug | 10,000 | 100 kW/m^2 | passive scalar thermal sensitivity, done as `v7_scalar_temperature_water_re10000_heatflux100kw` |
| B0 | FLiBe-like | 1,000 | 0 | serious-fluid low-Re debug |
| B1 | FLiBe-like | 10,000 | 0 | isothermal RANS baseline |
| B1a | FLiBe-like | 10,000 | 0 | isothermal RANS baseline, done as `v8_isothermal_rans_flibe_re10000_clean` |
| B2 | FLiBe-like | 10,000 | 100 kW/m^2 | passive thermal baseline, done as `v9_scalar_temperature_flibe_re10000_heatflux100kw` |
| B2r-flow | FLiBe-like | 10,000 | 0 | refined isothermal RANS repeat, done as `v10_isothermal_rans_flibe_re10000_refined` |
| B2r-thermal | FLiBe-like | 10,000 | 100 kW/m^2 | refined passive thermal repeat, done as `v12_scalar_temperature_flibe_re10000_heatflux100kw_refined_dt0005` |
| B2r2-flow | FLiBe-like | 10,000 | 0 | second thermal-focused RANS repeat, done as `v13_isothermal_rans_flibe_re10000_thermal_refined` |
| B2r2-thermal | FLiBe-like | 10,000 | 100 kW/m^2 | second thermal-focused passive repeat, done as `v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025` |
| B2r2-long | FLiBe-like | 10,000 | 100 kW/m^2 | passive long-run check, done as `v15_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_longrun_dt00025` |
| B3 | FLiBe-like | 10,000 | 500 kW/m^2 | high heat-flux sensitivity, done as `v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025` |
| C1 | FLiBe-like | 10,000 | 100 kW/m^2 | unsteady/LES attempt |

## Inlet Velocity Estimates

From `scripts/estimate_flow_params.py`, using the draft inlet size:

| Fluid | Re | Inlet Velocity | Mass Flow |
|---|---:|---:|---:|
| water debug | 1,000 | 0.0335 m/s | 0.0267 kg/s |
| water debug | 10,000 | 0.3348 m/s | 0.2670 kg/s |
| FLiBe-like | 1,000 | 0.1160 m/s | 0.1800 kg/s |
| FLiBe-like | 10,000 | 1.1598 m/s | 1.8000 kg/s |
| FLiBe-like | 50,000 | 5.7990 m/s | 9.0000 kg/s |

## Completed Run Ladder

1. A0: water, Re = 1,000, no heat. Done as `v3_isothermal_rans_water_re1000_corrected`.
2. A1: water, Re = 10,000, no heat. Done as `v5_isothermal_rans_water_re10000_clean`.
3. A2: water, Re = 10,000, small heat flux. Done as `v6_scalar_temperature_water_re10000_heatflux`.
4. A3: water, Re = 10,000, 100 kW/m^2. Done as `v7_scalar_temperature_water_re10000_heatflux100kw`.
5. B1a: FLiBe-like, Re = 10,000, no heat. Done as `v8_isothermal_rans_flibe_re10000_clean`.
6. B2: FLiBe-like, Re = 10,000, 100 kW/m^2. Done as `v9_scalar_temperature_flibe_re10000_heatflux100kw`.
7. B2r/B2r2: refined FLiBe-like RANS/passive mesh checks. Done as `v10`/`v12` and `v13`/`v14`.
8. B2r2-long: 100 kW/m^2 passive long-run check. Done as `v15_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_longrun_dt00025`.
9. B3: 500 kW/m^2 heat-flux sensitivity. Done as `v16_scalar_temperature_flibe_re10000_heatflux500kw_thermal_refined_dt00025`.

## Why This Order

Water debug makes solver mistakes easier to find.

FLiBe-like baseline gives the project fusion-blanket relevance.

The LES/unsteady case should not begin until B2 is stable, post-processed, and checked against a near-wall mesh/refinement sensitivity.

## Current Next Runs

| Priority | Case | Purpose |
|---:|---|---|
| 1 | advisor/report packaging | turn the completed v9/v12/v14/v15/v16 ladder into figures, assumptions, and limitations |
| 2 | one more targeted refinement | only if wall temperature becomes a central quantitative claim |
| 3 | C1 unsteady/LES attempt | capstone/stretch case after RANS/passive baseline is documented |

Completed reporting artifact:

- `docs/mesh_sensitivity_summary.md` presents `v9`/`v12`/`v14` as the first mesh-sensitivity result.
- `docs/passive_longrun_summary.md` presents `v14`/`v15` as the first time-extension result.
- `docs/heat_flux_sensitivity_summary.md` presents `v14`/`v16` as the first heat-load sensitivity result.
