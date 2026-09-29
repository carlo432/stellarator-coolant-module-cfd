# Stellarator First-Wall Coolant Module — CFD & Conjugate Heat Transfer

> **Status: ongoing.** This is an active individual research project. Results and documents are
> working drafts and will change.

A thermal-hydraulic feasibility study of a simplified, stellarator-relevant FLiBe coolant module, built
in **OpenFOAM** with **Gmsh** geometry. The project climbs a fidelity ladder:
steady RANS → passive-scalar heat transfer → unsteady RANS → wall-modelled LES →
wall-resolved conjugate heat transfer (CHT). It then couples in plasma-facing heat loads and
OpenMC volumetric heating.

It is a method-application study. It does not model a named reactor: the module is a controlled
geometry that isolates how coolant flow structure governs wall temperature.

## Central finding

> **Wall-temperature risk follows near-wall flushing.** Slow, broad near-wall flow creates
> hotspot-prone regions. Faster flushing reduces wall superheat, but hydraulic pumping cost rises steeply.

- **Geometry acts through flushing.** Outlet-offset and backward-facing-step geometries create slow
  near-wall flow and are hotspot-prone. A 90° bend is *not* automatically hot: its peak (59.9 K) is
  below both the straight baseline (63.5 K) and the outlet-offset case (77.0 K).
- **Wall-resolved CHT sets the physical temperature scale:** 27–28 K interface superheat at
  100 kW/m². Two meshes agree to within 3.4 %, and Dittus–Boelter / Gnielinski agree to within 6–10 %.
- **Pumping is the design knob.** Doubling Re from 10,000 to 20,000 drops interface ΔT from 28.2 K to
  17.0 K, while hydraulic pumping power rises ~7× (1.6 W → 11.5 W). Across the sweep ΔT ∝ U^−0.70 and
  pumping ∝ U^2.8.
- **Property realism matters.** Temperature-dependent FLiBe μ(T), k(T) lower hotspot superheat by
  ~18 % compared with constant properties.
- **Single-cell peak temperature is not a trustworthy metric.** Mesh studies show it is non-monotonic
  under refinement, so peak values are reported as same-mesh diagnostics only.

<p align="center">
  <img src="results/report_figures/35_central_finding_onepager.png" width="85%" alt="Central finding">
</p>

## Fidelity ladder

| Stage | Cases | Outcome |
|---|---|---|
| Geometry and mesh | `baseline_module` v1–v6 `.geo` | Gmsh module; best mesh 73k cells with 3,216 heated-wall faces |
| Steady RANS | v8–v13 | FLiBe-like Re = 10,000 baseline, Δp ≈ 1.77 kPa |
| Passive heat transfer | v9–v16 | 100 / 500 kW/m² heat-flux sweep and 3-level mesh sensitivity |
| Unsteady RANS | v17–v18, v25 | stays close to steady RANS |
| Geometry variants | v22–v31 | outlet-offset, 90° bend, backward-facing step |
| LES diagnostics | v26–v27 | coarse wall-modelled WALE LES and thermal LES (illustrative, not validated) |
| Conjugate heat transfer | CHT bridge, P4 geometry table, velocity sweeps | wall-resolved k-ω SST; correlation-validated |
| Channel LES method | periodic channel at Re_τ ≈ 171, Pr = 5 | qualified method check against DNS statistics |
| Coupled loads | OpenMC FLiBe heating, plasma heat-flux maps | mapped volumetric-source CHT |
| ARC precedent | 2-D ARC liquid-immersion blanket | OpenFOAM rebuild of published ARC blanket CFD (Ferrero; Leffler) |

The integrated write-up is [docs/final_report_draft_v0.md](docs/final_report_draft_v0.md).

## Repository layout

```
cases/
  baseline_module/     Gmsh .geo geometry and OpenFOAM case setups (v2–v31)
  branch_d/            CHT, LES and coupled-load cases, with their generator and
                       analysis scripts (see below)
  arc_blanket_2d*/     2-D ARC liquid-immersion blanket rebuilds
  hpc1_channel_*/      periodic-channel LES mesh/decomposition preflight
scripts/               case builders, post-processing, GCI, correlation checks, plotting,
                       and SLURM templates for HPC runs
results/
  report_figures/      numbered report figures
  */                   per-study CSV data and plots
docs/                  report draft, central finding, CHT validation, design notes
references/            bibliography (references.bib)
```

Case directories hold setup files only (`system/`, `constant/`, `0/`, geometry sources). Meshes, time
directories, solver logs and files over 1 MB are excluded for size; regenerate them from the `.geo`
files with Gmsh and `gmshToFoam`, or with the `gen_*.py` / `make_*` scripts in `cases/branch_d/`.

`cases/branch_d/` groups the higher-fidelity work:

| Cases | Study |
|---|---|
| `v22`–`v31` | outlet-offset, 90° bend and backward-facing-step geometries (RANS, URANS, LES) |
| `cht_*`, `cht_wr_U*` | conjugate heat transfer: mesh study, wall-resolved bridge, velocity sweep |
| `c2_*`, `c7_*`, `c13_*` | curved first-wall slice CHT, curved velocity sweep, curvature sweep |
| `c3_*`, `openmc_*` | OpenMC volumetric heating mapped into CHT |
| `geom_*` | plasma-load pipeline: shaped wall loads, phase and sensitivity sweeps, tritium, MHD drag |
| `a1_*`, `a2_*`, `a3_*`, `c1_*` | 3-D OpenMC torus sector, coil/SOL surrogate, variable-property CHT, wall-load inverse design |
| `les_channel_*`, `c9_thermal_gateclose`, `c14_*` | periodic-channel LES: momentum benchmarks, Pr = 5 thermal certification, Pr = 14.4 transfer check |

The Pr = 5 thermal-LES check compares against Kawamura et al. (1998) DNS profiles digitized from the
paper's figures (`kawamura_pr5_dns_digitized/`), with an estimated reading accuracy of about ±1 Θ⁺ unit.
It is not a comparison against tabulated DNS data.

## Requirements

- OpenFOAM v2512 (ESI)
- Gmsh 4.8
- Python 3.12 with numpy, scipy, matplotlib, pandas
- OpenMC (coupled-heating studies only)

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
cd cases/baseline_module/openfoam_cases/v14_scalar_temperature_flibe_re10000_heatflux100kw_thermal_refined_dt00025
# generate the mesh from cases/baseline_module/geometry, then run the solver named in system/controlDict
```

## Limitations

- The module is simplified and uses FLiBe-like properties. It is not a blanket-scale design.
- LES results are coarse diagnostics. Production wall-resolved LES is planned on HPC.
- MHD, tritium transport and structural response are out of scope or only scoped.
- Pumping power is hydraulic `Δp·Q` only.

## Acknowledgements

Geometry and precedent come from published ARC blanket CFD work
(Ferrero; Leffler); full citations are in `references/references.bib`. AI coding assistants were used
to help with implementation and to independently check results.

## License

Code is released under the MIT license (see [LICENSE](LICENSE)). Documents and figures are © the author.
