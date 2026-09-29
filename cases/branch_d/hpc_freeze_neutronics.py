#!/usr/bin/env python3
"""STEP 5 of the HPC streamlining plan: FREEZE the OpenMC heating field into a reusable OpenFOAM source.

OpenMC<->CFD is ONE-WAY (neutron/gamma heating -> flow) and the heating field is geometry-fixed, so run OpenMC
ONCE, then attach the frozen q'''(depth) as an fvOptions source to EVERY CFD tier -- never re-run neutronics per
CFD iteration. (C3 proved the path: OPENMC_CROSS_SECTIONS=/root/openmc_data; exe /root/openmc-src/build/bin/openmc.)

Reads an OpenMC heating-profile CSV in the C3 format (columns include x_mid_m, normalized_qvol_W_m3), fits
q'''(depth), and writes an fvOptions block:
  - mode=uniform : single scalarSemiImplicitSource, volumeMode absolute = total power (good when ~flat).
  - mode=mapped  : scalarSemiImplicitSource exprField, volumeMode specific = q'''(x) polynomial in depth.

Usage: hpc_freeze_neutronics.py <heating_csv> [--mode uniform|mapped] [--field h] [--coord x]
                                 [--total-power W] [--out fvOptions.frozen]
"""
import argparse, csv, sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("csv")
ap.add_argument("--mode", choices=["uniform", "mapped"], default="auto")
ap.add_argument("--field", default="h")
ap.add_argument("--coord", default="x", help="depth coordinate component: x|y|z")
ap.add_argument("--total-power", type=float, default=None, help="override total power [W]; else integrated")
ap.add_argument("--zones", type=int, default=8, help="mapped mode: number of piecewise-uniform depth zones")
ap.add_argument("--out", default=None)
args = ap.parse_args()

rows = list(csv.DictReader(open(args.csv)))
xs = [float(r["x_mid_m"]) for r in rows]
q = [float(r["normalized_qvol_W_m3"]) for r in rows]
dx = [float(r["x1_m"]) - float(r["x0_m"]) for r in rows]
qmean = sum(q) / len(q)
qpeak = max(q)
flatness = qpeak / qmean
# total power per unit cross-section area (W/m2) = integral q''' dx ; full power needs the wetted area
power_per_area = sum(qi * dxi for qi, dxi in zip(q, dx))

mode = args.mode
if mode == "auto":
    mode = "uniform" if flatness < 1.05 else "mapped"

print(f"=== freeze neutronics: {args.csv} ===")
print(f"  bins={len(rows)}  q'''_mean={qmean:.4g} W/m3  peak/mean={flatness:.4f}  -> mode={mode}")
print(f"  integral q''' dx = {power_per_area:.4g} W/m2 (multiply by wetted area for total power)")

pos = f"pos().{args.coord}()"
if mode == "uniform":
    tp = args.total_power if args.total_power is not None else None
    val = tp if tp is not None else f"<TOTAL_POWER_W>  # = {qmean:.6g} * fluid_volume; or pass --total-power"
    block = f"""frozenNeutronicsHeating
{{
    type            scalarSemiImplicitSource;
    active          true;
    selectionMode   all;
    volumeMode      absolute;            // value is TOTAL power [W]
    sources {{ {args.field} ( {val} 0 ); }}
}}"""
else:
    # PIECEWISE-UNIFORM DEPTH ZONES -- the exprField fvOption does NOT exist in OF2512 and scalarCodedSource
    # cannot dlopen as root (both proven by feature_smoke.py). Approximate q'''(depth) as a staircase of
    # cellZones, each driven by the VERIFIED scalarSemiImplicitSource. Emits a topoSetDict to build the zones.
    import numpy as np
    M = min(len(rows), args.zones)
    edges = np.linspace(min(float(r["x0_m"]) for r in rows), max(float(r["x1_m"]) for r in rows), M + 1)
    ax = {"x": 0, "y": 1, "z": 2}[args.coord]
    bands = []
    for b in range(M):
        lo, hi = edges[b], edges[b + 1]
        inb = [qi for qi, xi in zip(q, xs) if lo - 1e-12 <= xi < hi + 1e-12]
        bands.append((lo, hi, sum(inb) / len(inb) if inb else 0.0))
    def box(lo, hi):
        mn, mx = [-1e6, -1e6, -1e6], [1e6, 1e6, 1e6]; mn[ax], mx[ax] = lo, hi
        return f"({mn[0]:g} {mn[1]:g} {mn[2]:g}) ({mx[0]:g} {mx[1]:g} {mx[2]:g})"
    topo = "FoamFile { version 2.0; format ascii; class dictionary; object topoSetDict; }\nactions\n(\n"
    for b, (lo, hi, qb) in enumerate(bands):
        topo += (f"    {{ name qzone{b}; type cellSet; action new; source boxToCell; sourceInfo {{ box {box(lo,hi)}; }} }}\n"
                 f"    {{ name qzone{b}; type cellZoneSet; action new; source setToCellZone; sourceInfo {{ set qzone{b}; }} }}\n")
    topo += ");\n"
    block = "\n".join(
        f"""qzone{b}_heat
{{
    type            scalarSemiImplicitSource;
    active          true;
    selectionMode   cellZone;
    cellZone        qzone{b};
    volumeMode      specific;            // q''' [W/m3] frozen from OpenMC, band-averaged
    sources {{ {args.field} ( {qb:.6g} 0 ); }}
}}""" for b, (lo, hi, qb) in enumerate(bands))

out = args.out or (Path(args.csv).stem + ".fvOptions.frozen")
Path(out).write_text("FoamFile { version 2.0; format ascii; class dictionary; object fvOptions; }\n" + block + "\n")
print(f"  wrote {out}  (mode={mode}, field={args.field})")
if mode == "mapped":
    topo_out = Path(out).with_suffix(".topoSetDict")
    Path(topo_out).write_text(topo)
    print(f"  wrote {topo_out}  ({M} depth zones) -> run `topoSet -dict {topo_out.name}` BEFORE the solver")
print("  -> attach to a CFD tier's fluid-region fvOptions; run OpenMC ONCE, reuse across all tiers.")
