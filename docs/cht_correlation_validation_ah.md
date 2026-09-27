# CHT Correlation Validation

Date: 2026-06-20

Roadmap item: AH

Purpose: validate the wall-resolved CHT wall-temperature scale against independent turbulent-duct heat-transfer correlations.

## Bottom Line

The wall-resolved CHT bridge result is now the strongest wall-temperature result in the package:

> At `100 kW/m^2`, the wall-resolved kOmegaSST CHT bridge case gives interface superheat `27-28 K`. Two wall-resolved meshes differ by `3.4%`, solid conduction closes exactly, and the result agrees with Dittus-Boelter/Gnielinski film-temperature estimates within about `6-10%`.

This makes the passive-scalar peak wall `dT` a retracted/naive hotspot metric, not a trusted wall-temperature value.

## Inputs

| Quantity | Value |
|---|---:|
| Heat flux `q''` | `100,000 W/m^2` |
| Reynolds number | `10,000` |
| Prandtl number | `14.4` |
| Hydraulic diameter `Dh` | `0.02667 m` |
| Fluid thermal conductivity `k` | `1.0 W/m/K` |
| Gnielinski friction factor input | `0.0315` |
| Bulk rise over heated channel | about `0.37 K` |

Because the bulk rise is small, interface-minus-inlet is a clean approximation to film `dT` for this bridge check.

## Validation Table

Generated source:

- `results/branch_d/cht_correlation_validation.csv`
- `results/branch_d/cht_correlation_validation.png`
- `results/report_figures/31_cht_correlation_validation.png`

| Metric | `dT` | `Nu` | Interpretation |
|---|---:|---:|---|
| wall-resolved CHT, 50k | 28.16 K | 94.7 | within `9.5%` of Gnielinski film `dT` |
| wall-resolved CHT, 98k | 27.21 K | 98.0 | within `5.8%` of Gnielinski film `dT` |
| Dittus-Boelter | 25.17 K | 105.9 | independent turbulent-duct film scale |
| Gnielinski | 25.71 K | 103.7 | independent turbulent-duct film scale |
| passive peak, medium | 77.00 K | 34.6 | about `3x` above the correlation film scale |
| passive peak, fine | 105.38 K | 25.3 | about `4x` above the correlation film scale |

## Formulae Used

CHT-derived Nusselt number:

```text
Nu_CHT = (q'' / dT_CHT) * Dh / k
```

Dittus-Boelter:

```text
Nu = 0.023 * Re^0.8 * Pr^0.4
```

Gnielinski:

```text
Nu = [(f/8) * (Re - 1000) * Pr] /
     [1 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1)]
```

Film temperature difference from a correlation:

```text
dT = q'' * Dh / (Nu * k)
```

## Interpretation

The validation closes the loop that the passive-scalar model left open:

- passive scalar peak `dT`: non-monotonic and `3-4x` above correlation scale
- k-epsilon CHT: physical but wall-function scattered because the tested meshes sit in the buffer layer
- wall-resolved CHT: `27-28 K`, two meshes `3.4%` apart, correlation-consistent

Report-safe wording:

> The trusted wall-temperature scale in the current package starts from the wall-resolved CHT bridge result, `27-28 K` interface superheat at `100 kW/m^2`. It is grid-confirmed by two wall-resolved meshes and independently consistent with Dittus-Boelter/Gnielinski film-temperature estimates. The later P4 straight/BFS/outlet-offset CHT table extends this into a simplified-module geometry mechanism check; the passive Branch D four-geometry peaks remain diagnostic context, not physical wall-temperature predictions.

Citation note: final report should cite a heat-transfer textbook or original sources for Dittus-Boelter and Gnielinski before submission.
