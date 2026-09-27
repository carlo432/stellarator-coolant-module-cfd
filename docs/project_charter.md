# Project Charter

## Working Name

Revengeance

## Real Goal

Build a simplified stellarator-relevant 3D cooling/blanket module with heat transfer, run a RANS baseline, and attempt one LES or unsteady higher-fidelity case.

## Undersold Initial Pitch

Feasibility study of coolant-flow recirculation and heat-removal behavior in a simplified stellarator-relevant module.

This wording is intentionally conservative. It does not promise a full stellarator blanket, full MHD, tritium transport, or neutronics.

## What This Is Building On

Christopher Leffler's ARC blanket project:

- simplified 2D ARC blanket geometry
- Nek5000 LES
- comparison against existing RANS results
- focus on recirculation, split flow, and separation
- temperature and 3D left as future work

Our next-step-above target:

- simplified stellarator-relevant 3D geometry
- baseline RANS
- heat transfer
- one LES/unsteady case if feasible
- focus on recirculation, pressure drop, and hotspot risk

## Minimum Successful Project

The project succeeds if it produces:

1. A defensible simplified 3D coolant geometry.
2. A mesh with documented quality.
3. A working RANS baseline.
4. A heat-transfer case with an imposed heat flux or volumetric heating.
5. Plots showing velocity, pressure, temperature, and recirculation/stagnation regions.
6. A written comparison explaining what the model can and cannot claim.

## Strong Successful Project

The project is strong if it also includes:

1. A coarse LES or unsteady simulation.
2. RANS vs LES/unsteady comparison.
3. Line profiles matching the style of Leffler's post-processing.
4. Sensitivity study over inlet speed, heat flux, or mesh resolution.
5. A clear proposal for adding MHD/neutronics later.

## Out Of Scope For Senior-Year Version

- Full stellarator blanket.
- Full device CAD.
- Coupled neutronics-CFD-MHD-tritium.
- Validated reactor design claims.
- Tritium breeding ratio.
- Material stress/lifetime analysis.

## Key Research Question

Can a simplified stellarator-relevant cooling module show recirculation or low-flow regions that meaningfully affect heat removal, and does a higher-fidelity/unsteady simulation predict those regions differently than a RANS baseline?

## First Technical Decision Needed

Choose the geometry basis:

- W7-X-inspired module: easier public context, but W7-X is not a power reactor blanket.
- ARIES-CS-like module: more reactor-relevant and now backed by Ku compact-reactor, Ku physics-design, and Steiner divertor-summary source notes, but older and harder to simplify.
- Helios-like planar-coil stellarator module: modern and reactor-relevant, but still preconceptual.
- Fully idealized stellarator-relevant channel: easiest to execute, weakest connection to a named device.

Current recommended choice: idealized 3D module inspired by stellarator blanket/divertor constraints, with a short literature justification and explicit caveats that the module is not ARIES-CS CAD.
