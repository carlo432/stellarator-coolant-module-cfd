# FLiBe Thermal-Margin Chapter

Date: 2026-06-20

Roadmap item: AI

Purpose: convert the wall-resolved CHT bridge result into an operating-window statement for the current FLiBe-like coolant model.

## Bottom Line

At `100 kW/m^2`, the wall-resolved CHT bridge result gives a coolant-facing interface temperature of about `827-828 K` for an `800 K` inlet. Adding the exact `q*t/kappa` conduction drop through the `4 mm`, `kappa = 30 W/m/K` solid gives a heated solid-base temperature of about `841 K`. That is well above the FLiBe melt/liquidus guide near `732 K`, far below the broad reported liquid upper range near `1703 K`, and below the `873 K` upper end of the RAFM/Pint low-temperature structural-material window by about `31 K`.

At `500 kW/m^2`, The completed `cht_wallres_500` wall-resolved CHT run gives a coolant-facing interface temperature of `940.80 K` and a heated solid-base temperature of `1007.47 K`. That is not a safe-design claim. It is a computed warning case: the interface lies inside Pint's cited `600-700 deg C` (`873-973 K`) FLiBe blanket peak-operating band, but the solid base exceeds that band and is well above the nominal RAFM upper-temperature guide.

## Source Anchors

| Quantity | Value used here | Source note |
|---|---:|---|
| FLiBe melt/liquidus guide | `732 K` | `references/source_notes/attarian_flibe_thermophysical_properties.md`; `references/source_notes/pint_materials_assessment_flibe_fusion_blankets.md` |
| Reported broad FLiBe liquid range | about `732-1703 K` | `references/source_notes/attarian_flibe_thermophysical_properties.md` |
| FLiBe blanket peak operating range discussed by Pint | `600-700 deg C` = `873-973 K` | `references/source_notes/pint_materials_assessment_flibe_fusion_blankets.md` |
| RAFM steel operating-temperature guide | about `550-600 deg C` = `823-873 K` | `references/source_notes/pint_materials_assessment_flibe_fusion_blankets.md` |
| Current inlet temperature | `800 K` | project baseline |
| Wall-resolved CHT interface superheat | `27.21-28.16 K` at `100 kW/m^2`; `140.80 K` at `500 kW/m^2` | `docs/cht_correlation_validation_ah.md`; `docs/wall_resolved_500kw_cht_aj.md` |
| Heated-channel bulk rise | `0.37 K` at `100 kW/m^2` | CHT energy balance, `Q/(mdot*cp)` |
| Solid base-interface drop | `13.333 K` at `100 kW/m^2` | exact `q*t/kappa = 100000*0.004/30` |

## Thermal-Margin Table

Generated source:

- `results/branch_d/flibe_thermal_margin_ai.csv`
- `results/branch_d/flibe_thermal_margin_ai.png`
- `results/report_figures/32_flibe_thermal_margin.png`

| Case | Status | Interface temperature | Solid-base temperature | Margin statement |
|---|---|---:|---:|---|
| `100 kW/m^2` wall-resolved CHT bridge | run | `827.21-828.16 K` | `840.55-841.49 K` | above melt by about `95 K` at the interface and `109 K` at the base; solid base remains about `31 K` below `873 K` |
| `500 kW/m^2` wall-resolved CHT bridge | run | `940.80 K` | `1007.47 K` | interface enters the Pint blanket peak-operating band; solid base exceeds `973 K` by about `34 K` |

## Calculation Convention

Interface absolute temperature is taken directly from the wall-resolved CHT patch averages:

```text
T_interface = areaAverage(T) on heater_to_bottomWater
```

The heated solid-base temperature is taken from the matching solid patch average:

```text
T_base = areaAverage(T) on heatedBase
```

The exact conduction drop remains a check:

```text
T_base - T_interface ~= q'' * t_solid / k_solid
```

For the `500 kW/m^2` row, the completed CHT run confirms the expected linear scaling from the `100 kW/m^2` bridge:

```text
140.80 K / 28.16 K ~= 5.00
```

## Advisor-Safe Interpretation

The current `100 kW/m^2`, `800 K` inlet case is thermally plausible as a first FLiBe-like CHT bridge: it stays liquid, the trusted interface superheat is only `27-28 K`, and the heated solid base remains below the `873 K` upper end of the RAFM/Pint low-temperature material window. The margin is not huge, though. The computed `500 kW/m^2` wall-resolved CHT case puts the solid base at `1007.47 K`, above both the nominal RAFM upper guide and Pint's cited `973 K` upper end of the FLiBe blanket peak-operating band.

## Caveats

- This AI margin statement is still anchored to the straight-channel CHT bridge; use the later P4 table for the straight/BFS/outlet-offset geometry mechanism check.
- The current model uses constant FLiBe-like properties. Attarian and Sohal both support the limitation that FLiBe properties vary with temperature.
- The material-temperature discussion is a screening statement only. It does not solve corrosion, irradiation, redox chemistry, activation, tritium, or structural stress.
- The broad `732-1703 K` FLiBe liquid range is not the design window. Structural-material and blanket-material constraints are tighter in the present project.
