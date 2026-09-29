#!/usr/bin/env python3
"""C7: fit the curved-slice flushing law (interface wall dT vs inlet speed) and compare the exponent to the
straight-channel dT ~ U^-0.70. Reads the postprocessor summary JSONs written by the sweep driver."""
import json, math
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORK = Path(__file__).resolve().parent
FIGS = WORK / "figs"
CASES = [
    ("c7_curved_sweep_u058", 0.58), ("c7_curved_sweep_u087", 0.87),
    ("c7_curved_sweep_u116", 1.16), ("c7_curved_sweep_u174", 1.74),
    ("c7_curved_sweep_u232", 2.32),
]
NU = 3.09278e-6
DH = 0.0208696  # C2 curved hydraulic diameter (from C2 closure doc)

rows = []
for case, u in CASES:
    j = FIGS / f"{case}_summary.json"
    if not j.exists():
        print(f"MISSING {j}"); continue
    d = json.loads(j.read_text())
    tin = d["tin_K"]
    iface = d["patch_averages"]["heater_to_bottomWater"]["area_average"]
    base = d["patch_averages"]["heatedBase"]["area_average"]
    dT_iface = iface - tin
    dT_base = base - tin
    Re = u * DH / NU
    rows.append((u, Re, dT_iface, dT_base, d["solver"]["execution_time_s"],
                 d["solver"].get("fluid_h_final_residual_final_step")))
    print(f"u_in={u:.2f}  Re={Re:7.0f}  interface dT={dT_iface:6.3f} K  base dT={dT_base:6.3f} K")

if len(rows) >= 2:
    us = [r[0] for r in rows]; dTs = [r[2] for r in rows]
    lx = [math.log(u) for u in us]; ly = [math.log(t) for t in dTs]
    n = len(lx); sx = sum(lx); sy = sum(ly); sxx = sum(x*x for x in lx); sxy = sum(x*y for x, y in zip(lx, ly))
    slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)       # exponent
    inter = (sy - slope*sx) / n
    ybar = sy/n; sstot = sum((y-ybar)**2 for y in ly)
    ssres = sum((y-(inter+slope*x))**2 for x, y in zip(lx, ly))
    r2 = 1 - ssres/sstot if sstot else float("nan")
    print("=" * 64)
    print(f"  CURVED flushing law:  interface dT ~ U^{slope:.3f}   (R^2={r2:.4f})")
    print(f"  STRAIGHT reference :  interface dT ~ U^-0.70")
    print(f"  -> curved exponent {'STEEPER' if slope < -0.70 else 'SHALLOWER'} than straight by {abs(slope+0.70):.3f}")
    print("=" * 64)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(us, dTs, "o", ms=8, color="tab:green", label="curved CHT (wall-resolved)")
    xf = [min(us)*0.95, max(us)*1.05]
    ax.loglog(xf, [math.exp(inter)*x**slope for x in xf], "-", color="tab:green",
              label=f"curved fit  dT~U^{slope:.2f} (R²={r2:.3f})")
    a070 = dTs[2] / us[2]**(-0.70)
    ax.loglog(xf, [a070*x**(-0.70) for x in xf], "--", color="tab:red", label=r"straight $dT\sim U^{-0.70}$")
    ax.set_xlabel("inlet speed U [m/s]"); ax.set_ylabel("interface wall ΔT [K]")
    ax.set_title("C7: flushing law on curvature vs straight channel")
    ax.grid(True, which="both", alpha=0.3); ax.legend()
    out = FIGS / "c7_curved_flushing_law.png"; fig.tight_layout(); fig.savefig(out, dpi=140)
    print("wrote", out)

    # csv
    csv = FIGS / "c7_curved_flushing_law.csv"
    with open(csv, "w") as f:
        f.write("u_in,Re,interface_dT_K,base_dT_K,exec_s\n")
        for r in rows: f.write(f"{r[0]},{r[1]:.0f},{r[2]:.4f},{r[3]:.4f},{r[4]}\n")
    print("wrote", csv)
