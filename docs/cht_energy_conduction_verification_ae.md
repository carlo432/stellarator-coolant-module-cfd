# CHT Energy-Balance And Solid-Conduction Verification

Date: 2026-06-20

Roadmap item: AE

Purpose: document the basic physical closure checks for the CHT bridge cases.

## Bottom Line

The CHT bridge cases pass the most important local conservation check: the solid temperature drop from heated base to fluid interface matches the analytic conduction drop

```text
q'' * t / kappa = 100000 * 0.004 / 30 = 13.333 K
```

to roundoff for every retained CHT mesh. In `cht_v1`, the fluid-side and solid-side interface patch averages also match exactly at `843.544745 K`, confirming the coupled-interface temperature continuity for that k-epsilon medium-mesh case.

Important interpretation: `cht_v1` is not the converged wall-temperature result. It is the superseded k-epsilon medium-mesh/wall-function case. The current trusted wall-temperature scale is the wall-resolved kOmegaSST CHT bridge: interface temperature about `827-828 K`, or `27-28 K` above the `800 K` inlet. The `cht_v1` value remains useful for checking interface continuity and solid conduction, not for reporting the final wall temperature.

The available `cht_v1` outlet value is an area-average outlet temperature rise of `0.427 K`. The expected bulk energy-balance rise is `Q/(mdot*cp) = 1600/(1.80*2400) = 0.370 K`. This difference should be treated the same way as the earlier outlet-temperature caveat: area-average outlet temperature is not a mass-flow-weighted bulk outlet metric. It is close in scale, but a true CHT bulk-energy closure should use a mass-flow-weighted outlet post-process.

## Generated Source

- `results/branch_d/cht_energy_conduction_verification_ae.csv`

## Verification Table

| Case | Model | Interface dT | Base-interface drop | Expected drop | Verification |
|---|---|---:|---:|---:|---|
| `cht_coarse` | k-epsilon CHT | `23.94 K` | `13.333351 K` | `13.333333 K` | solid conduction closes |
| `cht_medium` / `cht_v1` | k-epsilon CHT | `43.54 K` | `13.333326 K` | `13.333333 K` | solid conduction closes; interface matches both sides |
| `cht_fine` | k-epsilon CHT | `28.61 K` | `13.333324 K` | `13.333333 K` | solid conduction closes |
| `cht_wallres` | wall-resolved kOmegaSST CHT | `28.16 K` | `13.333344 K` | `13.333333 K` | solid conduction closes |
| `cht_wallres_fine` | wall-resolved kOmegaSST CHT | `27.21 K` | `13.333312 K` | `13.333333 K` | solid conduction closes |

## Energy-Scale Check

The imposed heat input is:

```text
Q = q'' * A = 100000 W/m^2 * 0.016 m^2 = 1600 W
```

Using the CHT bridge mass flow and heat capacity:

```text
mdot = 1.80 kg/s
cp = 2400 J/kg/K
dT_bulk = Q / (mdot * cp) = 0.370 K
```

The available `cht_v1` outlet patch average is `800.427 K`, or `0.427 K` above the `800 K` inlet. That is the right small-temperature-rise scale, but it should not be called exact energy closure because the retained output is an area average. For report wording, use:

> The CHT bridge heat input implies a bulk coolant rise of only about `0.37 K`, so the wall/interface superheat is dominated by the local film resistance, not by bulk coolant warming. The retained outlet patch average is about `0.43 K`, consistent in scale but not a mass-flow-weighted bulk metric.

## What This Verifies

- The imposed solid-base gradient corresponds to the intended `100 kW/m^2` heat flux.
- The `4 mm`, `kappa = 30 W/m/K` solid slab produces the exact expected `13.333 K` base-to-interface drop.
- The coupled interface is continuous in `cht_v1`.
- The wall-resolved interface temperature scale, `27-28 K`, is not a numerical artifact of a broken solid conduction setup.

## What This Does Not Verify

- It does not make the k-epsilon wall-function CHT grid sequence converged.
- It does not replace the AH Dittus-Boelter/Gnielinski validation.
- By itself, it does not prove the later P4 geometry-level CHT mechanism table.
- It does not provide a mass-flow-weighted CHT outlet bulk-temperature closure yet.
