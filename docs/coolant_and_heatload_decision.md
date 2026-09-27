# Coolant And Heat-Load Decision

## Recommendation

Use an imposed heated wall and FLiBe-like molten salt as the serious case.

## Why Heated Wall First

A wall heat flux is easier to implement and explain than a volumetric neutron heating profile.

It also maps naturally to the question:

> Does coolant recirculation near a heated plasma-facing boundary create local heat-removal problems?

## Why Not Volumetric Heating First

Volumetric heating is more blanket-realistic, but it requires assumptions about neutron deposition and material regions. That pulls the project toward neutronics before the CFD baseline is stable.

Use volumetric heating as a future upgrade.

## Why FLiBe-Like Molten Salt

FLiBe-like molten salt keeps the project close to liquid blanket relevance while avoiding full liquid-metal MHD as the first physics stack.

Important caveat:

Molten salt still has material-property uncertainty, material-compatibility concerns, tritium chemistry, redox/corrosion issues, irradiation constraints, and some electromagnetic relevance. The 2024 ORNL FLiBe materials assessment supports FLiBe as a serious fusion breeder-coolant candidate, but it also makes clear that structural-material readiness is not solved. For this project, FLiBe-like coolant is a manageable first serious CFD baseline, not a blanket materials design.

The Sagara FFHR source adds a helical-reactor FLiBe blanket lineage and is now full-text local. Use it to support FLiBe precedent and a `0.1 MW/m^2` FFHR thermal/stress heat-flux scale, but not as validation of our geometry or wall temperatures.

The FLiBe high-Prandtl heat-transfer enhancement cluster adds a separate path for the later "switch once slices work" idea. Chiba 2005 and the wider Chiba/Shimizu/Watanabe/Kim cluster are now full-text local and support high-Pr molten-salt heat-transfer enhancement plus pressure-drop tradeoff language; extract the wider cluster only before claiming detailed packed-bed dimensions, correlations, or results.

## Debug Fluids

It is acceptable to debug with air or water.

That should be described as a numerical/debug phase, not the final physics case.

## FLiBe-Like Properties

These are first-pass constant properties for the current CFD baseline. They are now traceable to the ARC FLiBe property table in Sorbom et al. and cross-checked against independent molten-salt property sources, but they are still simplified and temperature-specific.

- density: 1940 kg/m^3
- dynamic viscosity: 0.006 Pa s
- thermal conductivity: 1.0 W/m K
- heat capacity: 2400 J/kg K
- inlet temperature: 800 K

Source note:

- `references/source_notes/sorbom_arc_flibe_blanket.md`
- `references/source_notes/pint_materials_assessment_flibe_fusion_blankets.md`
- `references/source_notes/tutwiler_twisted_tape_molten_salt_first_wall_2025.md`
- `references/source_notes/flibe_high_prandtl_heat_transfer_enhancement.md`
- `references/source_notes/sagara_ffhr_flibe_blanket.md`
- `references/source_notes/attarian_flibe_thermophysical_properties.md`
- `references/source_notes/sohal_inl_liquid_salt_properties.md`
- Local PDF: `references/pdfs/sorbom_arc_1409_3540.pdf`
- Local extracted text: `references/extracted_text/sorbom_arc_1409_3540.txt`

Important caveat:

Sorbom et al. list the FLiBe properties for liquid FLiBe at 950 K. Sohal's INL database gives viscosity and density correlations that place the current density and viscosity near the 900-950 K range, while Attarian et al. show that viscosity is higher near 800 K. The current OpenFOAM cases therefore use a constant-property FLiBe-like model with an 800 K inlet; this is appropriate for an early feasibility baseline but not a final temperature-dependent material model.

The ORNL materials assessment adds a second caveat: even when FLiBe is attractive as a breeder-coolant, compatibility with fusion-relevant low-activation structural materials, flowing-loop corrosion, redox control, tritium behavior, magnetic-field effects, and irradiation remain major open issues.

## First Heat Flux Values

- debug: 10 kW/m^2
- baseline: 100 kW/m^2
- stretch: 500 kW/m^2

Literature context:

- Sorbom ARC gives a first-wall/vacuum-vessel heat-load scale for the original ARC thermal analysis.
- Kuang ARC heat exhaust gives a much higher divertor heat-load scale, including MW/m^2-class divertor values.
- Mau et al.'s ARIES-CS divertor heat-load full text gives directly stellarator-specific target-plate context: a 10 MW/m^2 peak-load design goal and combined thermal/alpha peak loads around 5-18 MW/m^2 before further optimization.
- Sagara FFHR and Chiba 2005 both include `0.1 MW/m^2` heat-transfer/thermal-analysis scales, while Wang's ARIES-CS blanket scoping uses `0.5 MW/m^2` plasma heat flux.
- Steiner's short ARIES-CS final-report summary gives directly stellarator-specific divertor heat-load context, with a final-design thermal heat flux of 3.61 MW/m^2 under a 10 MW/m^2 prescribed limit.
- Lanahan et al. 2025 gives modern T-tube thermal-fluid/thermal-structural context for nonuniform and transient MW/m^2-class heat loads, but only metadata/abstract has been ingested locally.
- Tutwiler et al. 2025 gives modern molten-salt first-wall cooling context, describing `1-5 MW/m^2` incident heat-flux expectations for first-wall/divertor coolant channels and using LES to report Nusselt number and friction factor.
- Chiba 2005 and the wider FLiBe high-Prandtl cluster support future metric selection around heat-transfer enhancement, Nusselt/heat-transfer coefficients, friction/pressure-drop behavior, and sphere-packed channel concepts.
- Therefore, 100 and 500 kW/m^2 are best described as early wall-load sensitivity cases for the simplified module, not final reactor or divertor heat-flux claims.

Source notes:

- `references/source_notes/sorbom_arc_flibe_blanket.md`
- `references/source_notes/kuang_arc_heat_exhaust.md`
- `references/source_notes/tutwiler_twisted_tape_molten_salt_first_wall_2025.md`
- `references/source_notes/flibe_high_prandtl_heat_transfer_enhancement.md`
- `references/source_notes/mau_aries_cs_divertor_heat_load.md`
- `references/source_notes/steiner_aries_cs_divertor_final_report.md`
- `references/source_notes/lanahan_2025_ttube_thermal_fluid_structural.md`

## Decision To Confirm With Advisor

Ask whether the serious baseline should use:

- FLiBe-like molten salt
- helium
- water as an educational stand-in
- liquid metal with MHD deferred

Current recommendation: FLiBe-like molten salt.
