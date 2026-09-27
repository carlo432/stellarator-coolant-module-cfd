# C6 FLiBe Property Realism Status

Date: 2026-06-26

## Bottom Line

C6 is complete as a property-realism sensitivity note.  It does not add a
temperature-dependent OpenFOAM solve; it quantifies what the current
constant-property assumption means.

Report-safe claim:

> The current FLiBe-like kinematic viscosity corresponds to a high-temperature
> property point near `950 K`, not exact `800 K` FLiBe.  Over the local
> `~28 K` CHT wall-superheat scale, `nu` and `Pr` change by about `10%`, but
> across the wider `730-950 K` operating range they change by about a factor of
> three.  Forced convection dominates the tested baseline: `Gr/Re^2` is
> `6.571e-04` on the local CHT scale and about `0.013` on the later
> plasma-hotspot scale.  This is a scoped regime statement, not a general
> license to omit buoyancy.

## Correlation Result

The C6 script uses FLiBe property correlations:

- `rho = 2415.6 - 0.49072 T`
- `mu = 1.16e-4 exp(3755 / T)`
- `cp = 2386 J/kg/K`
- `k = 1.1 W/m/K`

The project LES value `nu = 3.09278e-6 m^2/s` implies:

| Quantity | Value |
|---|---:|
| implied reference temperature | `950.5 K` |
| density at reference | `1949.2 kg/m^3` |
| dynamic viscosity at reference | `0.00603 Pa s` |
| Pr at reference | `13.08` |

## Temperature Sensitivity

| T, K | rho, kg/m3 | mu, Pa s | nu, m2/s | Pr |
|---:|---:|---:|---:|---:|
| `730` | `2057.4` | `0.01988` | `9.6624e-06` | `43.12` |
| `800` | `2023.0` | `0.01267` | `6.2651e-06` | `27.49` |
| `850` | `1998.5` | `0.00962` | `4.8119e-06` | `20.86` |
| `900` | `1974.0` | `0.00752` | `3.8115e-06` | `16.32` |
| `950` | `1949.4` | `0.00604` | `3.0986e-06` | `13.10` |
| `960` | `1944.5` | `0.00580` | `2.9811e-06` | `12.57` |

Over the approximate `28 K` wall-superheat scale:

- `nu`: `3.093e-06 -> 2.782e-06`, or `-10.1%`
- `Pr`: `13.08 -> 11.68`, or `-10.7%`

Over `730 -> 950 K`:

- `nu` changes by `x3.12`
- `Pr` changes by `x3.29`

## Buoyancy Check

Using `Dh = 0.0208696 m`, `Re = 10000`, and `dT = 28 K`:

| Quantity | Value |
|---|---:|
| beta | `2.518e-04 1/K` |
| Gr | `6.571e+04` |
| Ri = Gr/Re^2 | `6.571e-04` |

Interpretation: `Ri = 6.571e-04`, so forced convection dominates this baseline.
The calculation bounds buoyancy at this operating point; it does not establish
that buoyancy can be omitted across the wider design space.

Follow-on 2026-07-09 plasma-hotspot audit:

| dT scale | Gr/Re^2 | margin to `0.1` flag |
|---|---:|---:|
| original C6 local CHT scale | `6.571e-04` | about `150x` |
| plasma hotspot, about `1223 K` superheat | `1.25e-02` | about `8x` |

The conclusion survives, but the margin is smaller for the plasma-shaped
hotspot.  Use the plasma-hotspot value when discussing the July pipeline load
envelope.

## Artifacts

- `cases/branch_d/c6_flibe_property_realism.py`
- `cases/branch_d/figs/c6_flibe_property_realism.png`

## Caveat

This closes C6 as a scoping calculation.

## A3 CFD Follow-On

A3 now supplies the first temperature-dependent OpenFOAM comparison requested by this caveat.
Matched plasma-load RANS-CHT cases with tabulated `mu(T)` and `k(T)` lower interface-average and
hotspot superheat by `10.9%` and `18.0%` relative to the legacy constants. Equivalently, the
legacy model biases those metrics high by `12.3%` and `21.9%`. See
`docs/a3_variable_property_cht_result.md`.

The same-speed/same-Re matrix remains useful as a broader operating-point study. A3 is one
same-speed representative channel comparison and also updates the conductivity reference level;
it is not a full property-validation campaign.
