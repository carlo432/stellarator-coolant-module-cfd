# Branch D Dimensionless Metrics Summary

Status: roadmap item F started and packaged; updated after the peak-dT GCI check.

This note translates the four-geometry Branch D flushing-speed result into diagnostic non-dimensional heat-transfer and loss metrics. After the outlet-offset GCI check showed non-monotonic single-cell peak dT, wall-average `Nu/St` should be treated as the preferred non-dimensional thermal metric. Peak-based `Nu` remains useful only as a same-mesh hotspot diagnostic.

## Artifacts

- Script: `scripts/plot_branch_d_dimensionless_metrics.py`
- Data: `results/branch_d/branch_d_dimensionless_metrics.csv`
- Figure: `results/branch_d/branch_d_dimensionless_metrics.png`
- Advisor copy: `results/report_figures/25_branch_d_dimensionless_metrics.png`

## Definitions Used

These are diagnostic passive-scalar metrics, not validated developed-channel correlations.

- FLiBe-like constants: `rho = 1940 kg/m^3`, `mu = 0.006 Pa s`, `k = 1.0 W/m/K`, `cp = 2400 J/kg/K`, `Pr = 14.4`
- Inlet hydraulic diameter: `Dh = 0.026667 m`
- Inlet speed: `U_inlet = 1.1598 m/s`
- Heat flux: `q'' = 100,000 W/m^2`
- Wall temperatures: VTP `heated_wall` surface cell data at the final exported time
- Preferred: `h_mean = q'' / mean_wall_deltaT`
- Preferred: `Nu_mean = h_mean * Dh / k`
- Preferred: `St_mean_local = h_mean / (rho * cp * U_nearwall)`
- Same-mesh hotspot diagnostic only: `h_peak = q'' / peak_wall_deltaT`
- Same-mesh hotspot diagnostic only: `Nu_peak = h_peak * Dh / k`
- `K_loss = 2 delta_p / (rho * U_inlet^2)`, where a documented or patch-averaged pressure drop is available
- `f_global = K_loss * Dh / L`, with `L = 0.40 m`, reported only as a global diagnostic

## Current Table

| Geometry | Case | Near-wall speed [m/s] | Mean wall dT [K] | Local Re_Dh | Nu_mean | St_mean_local x 1000 | Peak dT diagnostic [K] | K_loss |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 90-deg bend | v28/v29 | 0.860 | 43.51 | 7415 | 61.29 | 0.574 | 59.86 | 1.36 |
| straight | v13/v14 | 0.470 | 43.80 | 4052 | 60.89 | 1.043 | 63.52 | 1.36 |
| outlet-offset | v22/v23 | 0.390 | 44.16 | 3363 | 60.38 | 1.247 | 77.00 | 1.52 |
| backward-facing step | v30/v31 | 0.296 | 46.37 | 2552 | 57.51 | 1.565 | 77.44 | not reported |

## Interpretation

The wall-average metrics show that total/average wall response is much less geometry-discriminating than the single-cell peak. This matches the energy-balance story and the GCI check: wall average is a safer thermal quantity, while peak values identify same-mesh hotspot severity.

The outlet-offset and backward-facing-step cases still sit at the low-flushing/high-hotspot end in the peak diagnostic, but those peak values are mesh-sensitive. The defensible wording is: near-wall velocity fields support the flushing mechanism; peak temperature values are same-mesh hotspot indicators; wall-average `Nu/St` are the preferred non-dimensional thermal metrics until conjugate heat transfer or a grid-robust wall-temperature metric is available.

The BFS case remains the forward-prediction confirmation for the velocity/flushing trend, but its pressure-loss coefficient should not be quoted yet because no patch-averaged pressure-drop file is retained for v30. A naive cell-mean pressure proxy failed a sign check and is intentionally not used.

## Caveats

These metrics should be described as apparent, module-level diagnostics. They are not a validation against a canonical pipe, duct, bend, or BFS heat-transfer correlation. The thermal solve is passive scalar with constant properties, no buoyancy, no conjugate wall conduction, no MHD, no neutronics, no tritium transport, and no formal grid convergence for the geometry variants. The single-cell peak wall dT is explicitly mesh-sensitive after the outlet-offset GCI check; do not call the peak trend grid-converged.
