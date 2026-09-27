# Baseline Module Spec

## Purpose

Create the smallest 3D coolant module that can show separation, recirculation, pressure drop, and heat-removal behavior.

This is the first test geometry for the project.

## Geometry Concept

Expanded heated channel:

- Coolant enters through a narrow inlet channel.
- Flow expands into a wider heated chamber.
- Flow exits through a narrower outlet channel.
- The expansion should create separation/recirculation.
- One wall represents a plasma-facing heated surface.

## Draft Dimensions

All dimensions are placeholders for first testing.

- Inlet channel length: 0.10 m
- Inlet width: 0.02 m
- Chamber length: 0.20 m
- Chamber width: 0.08 m
- Outlet channel length: 0.10 m
- Outlet width: 0.02 m
- Depth: 0.04 m

Total streamwise length: 0.40 m

## Boundary Names

- `inlet`
- `outlet`
- `heated_wall`
- `adiabatic_walls`

## First Boundary Conditions

Isothermal RANS:

- Inlet: fixed velocity.
- Outlet: fixed pressure.
- Walls: no slip.

Thermal RANS:

- Inlet: fixed temperature.
- Heated wall: fixed heat flux.
- Other walls: adiabatic.

## First Fluid

Use water or air only for debug if needed.

For the serious case, choose one:

- FLiBe/molten salt for blanket relevance.
- Helium for gas-cooled simplicity.
- Liquid metal only if MHD becomes central.

## First Outputs

- Velocity magnitude.
- Pressure.
- Streamlines.
- Recirculation zone.
- Temperature.
- Max temperature.
- Pressure drop.

## Pass/Fail

The baseline geometry is useful if:

- it meshes cleanly
- the solver converges
- the expansion produces recirculation
- heat transfer creates a measurable temperature pattern
- post-processing can produce line profiles
