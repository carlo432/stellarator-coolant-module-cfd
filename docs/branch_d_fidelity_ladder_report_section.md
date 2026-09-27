# Branch D Fidelity Ladder Report Section

Date: 2026-06-20

Purpose: provide a paste-ready report section that connects the Branch D geometry comparison, RANS, passive thermal, heat-flux, URANS, first LES smoke-test, and thermal LES extension results.

## Section Draft

Branch D became the main comparison branch because it tests how geometry-driven near-wall flow structure affects local wall-temperature response while keeping the simplified module controlled. The branch now includes four geometries: the straight baseline, the outlet-offset recirculation case, and a 90-degree bend. The 90-degree bend is the key discriminating case because it contains a bend and local separation but does not produce the strongest hotspot.

Across the four tested geometries, the peak wall-temperature rise decreases as near-wall flushing speed increases: the 90-degree bend has near-wall mean speed about `0.86 m/s` and peak wall dT about `59.9 K`; the straight baseline has about `0.47 m/s` and `63.52 K`; the outlet-offset case has about `0.39 m/s` and `77.00 K`; and the backward-facing step has about `0.296 m/s` and `77.4 K`. The Pearson correlation is about `r = -0.85` for `n = 4`. This is not a validated correlation, but it gives a better working claim than "recirculation always makes hotspots." The BFS case confirms the slow-flow / hot-peak prediction.

The peak-temperature part of this trend is now explicitly same-mesh and mesh-sensitive. The outlet-offset GCI check was non-monotonic (`94.43 K`, `77.00 K`, `105.40 K` for coarse/medium/fine), so peak dT should not be called grid-converged. Wall-average metrics and velocity-side flushing evidence carry the quantitative weight.

The outlet-offset geometry keeps the same coolant model, Reynolds number, heated-wall area, heat flux, and wall discretization as the straight thermal-refined baseline, but shifts the outlet path so that the chamber flow turns more strongly near the heated wall.

The steady RANS result confirms that the outlet-offset geometry strengthens the low-flow/recirculation region. Relative to the straight baseline, the near-wall reverse-flow fraction increases from about `0.363` to `0.462`, the near-wall mean streamwise velocity drops from about `0.177 m/s` to `0.087 m/s`, and the pressure drop increases from about `1774 Pa` to about `1983 Pa`.

The passive thermal result shows that the effect is local rather than global. At `100 kW/m^2`, the wall-average temperature rise remains nearly unchanged (`43.80 K` straight versus `44.16 K` outlet-offset), but the peak wall-temperature rise increases from `63.52 K` to `77.00 K`. The peak also moves upstream into the recirculation region. This supports a local hotspot/distribution claim, not a validated blanket heat-removal prediction.

The 90-degree bend result prevents overgeneralization. Its peak wall dT is only `59.9 K`, below the straight and outlet-offset cases. The interpretation is that the bend's corner separation is weak/localized and reattaches; the wall-temperature peak is streamwise development, not broad recirculation-driven heating.

The heat-flux extension strengthens the same interpretation. At `500 kW/m^2`, the wall-average temperature rise remains nearly geometry-insensitive (`218.89 K` straight versus `220.83 K` outlet-offset), while the outlet-offset peak wall-temperature rise reaches `385.16 K` versus `316.61 K` for the straight case. The offset-straight peak gap therefore grows from about `+13.5 K` at `100 kW/m^2` to about `+68.6 K` at `500 kW/m^2`.

The unsteady RANS checks show that k-omega SST URANS does not create meaningful time dependence in this module. Straight v18 and bent v25 both remain close to their corresponding steady baselines. In the bent case, the recirculation structure is present, but probe fluctuations are very small, with residual velocity standard deviation about `0.0004 m/s`.

The first LES smoke test adds a model-hierarchy result. In v26, we ran a coarse, wall-modeled WALE LES on the bent geometry to `0.8 s` and judged the developed window after `t > 0.5 s` with detrending. The downstream recirculation probe, located near the v23 hotspot streamwise position, had residual velocity standard deviation about `0.054 m/s` and peak-to-peak fluctuation about `0.29 m/s`, far above the URANS residual. This suggests the hotspot-region recirculation can be dynamically active in LES even though URANS averages it steady.

The v27 thermal LES extension adds a first thermal-unsteadiness check in the same simplified module. The downstream/hotspot-region fluid probe about `1 mm` off the wall had detrended temperature fluctuation standard deviation about `0.234 K` and peak-to-peak variation about `1.49 K`. Resolved-field statistics also place the largest sampled mean resolved TKE in the hotspot-side recirculation band: about `0.0727 m^2/s^2` in the downstream hotspot band and about `0.0594 m^2/s^2` in the broader hotspot zone, with the inlet duct effectively zero.

The LES results must remain caveated. They are coarse wall-modeled LES cases on the RANS-scale mesh, not wall-resolved, not grid-converged, short in sample duration, and not validated. v26 is isothermal and does not update passive scalar wall-temperature values. v27 adds thermal fluctuations, but the promoted probe is near-wall fluid about `1 mm` off the wall, not the wall surface. The value of v26/v27 is that they qualitatively support the RANS-vs-LES motivation in the same simplified module.

## Recommended Figure

Use:

- `results/branch_d/branch_d_fidelity_ladder_panel.png`
- advisor copy: `results/report_figures/21_branch_d_fidelity_ladder_panel.png`
- near-wall flushing-speed trend: `results/report_figures/24_branch_d_nearwall_flushing_law.png`
- v27 mean-field comparison: `results/report_figures/22_v27_rans_vs_les_mean_velocity_midplane.png`
- v27 resolved-TKE and thermal-fluctuation view: `results/report_figures/23_v27_resolved_tke_temperature_fluctuation_midplane.png`

Caption:

> Branch D fidelity ladder for the outlet-offset geometry. The RANS velocity slice identifies the strengthened low-flow/recirculation region; the coarse LES-vs-URANS probe comparison shows recirculation-zone fluctuations that URANS averages steady; the wall-temperature comparison shows the local peak penalty at `100 kW/m^2`; and the heat-flux scaling panel shows the peak penalty growing from about `+13.5 K` at `100 kW/m^2` to about `+68.6 K` at `500 kW/m^2`. The LES panel is an illustrative smoke test only: coarse, wall-modeled, short, isothermal, not wall-resolved, not grid-converged, and not validated.
