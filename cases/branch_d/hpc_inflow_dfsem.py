#!/usr/bin/env python3
"""STEP 4 of the HPC streamlining plan: turn a developed PERIODIC channel (run once) into a synthetic-turbulence
INLET for the non-periodic geometry tiers (curved slice, stellarator segment) -- so they don't grow turbulence
over a long development length in the expensive domain.

Extracts U(y) mean profile + Reynolds-stress tensor R(y) from a developed channel's UMean/UPrime2Mean, and writes
OpenFOAM turbulentDFSEM (divergence-free synthetic-eddy method) inlet data:
  <out>/constant/boundaryData/<patch>/points
  <out>/constant/boundaryData/<patch>/0/{U,R,L}
plus the 0/U inlet BC snippet. One precursor feeds all geometry tiers; no plane-storage needed.

Usage: hpc_inflow_dfsem.py <developed_channel_case> [--time latestTime] [--patch inlet] [--delta auto]
                           [--out <dir>]
"""
import argparse, re, math, sys
from pathlib import Path
import numpy as np

def read_scalar(fn):
    s=open(fn).read(); m=re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)\s*;', s, re.S)
    return np.array([float(x) for x in m.group(1).split()])
def read_vector(fn):
    s=open(fn).read(); m=re.search(r'internalField\s+nonuniform\s+List<vector>\s*\d+\s*\((.*?)\)\s*;', s, re.S)
    return np.array([[float(v) for v in n.split()] for n in re.findall(r'\(([^()]*)\)', m.group(1))])
def read_symm(fn):
    s=open(fn).read(); m=re.search(r'internalField\s+nonuniform\s+List<symmTensor>\s*\d+\s*\((.*?)\)\s*;', s, re.S)
    return np.array([[float(v) for v in n.split()] for n in re.findall(r'\(([^()]*)\)', m.group(1))])

ap = argparse.ArgumentParser()
ap.add_argument("case")
ap.add_argument("--time", default=None)
ap.add_argument("--patch", default="inlet")
ap.add_argument("--delta", default="auto")
ap.add_argument("--out", default=None)
args = ap.parse_args()

case = Path(args.case)
t = args.time
if t is None:  # latest numeric time
    times = sorted([p.name for p in case.iterdir() if re.fullmatch(r"[0-9.]+", p.name)], key=float)
    t = times[-1]
Cy = read_scalar(case/t/"Cy"); Um = read_vector(case/t/"UMean"); R = read_symm(case/t/"UPrime2Mean")
delta = (Cy.max())/2 if args.delta == "auto" else float(args.delta)

ylevels = np.unique(np.round(Cy, 9))
def avg(arr, col=None):
    return np.array([(arr[np.isclose(Cy, y)] if col is None else arr[np.isclose(Cy, y), col]).mean() for y in ylevels])
yc = avg(Cy)
Uy = np.column_stack([avg(Um, 0), avg(Um, 1), avg(Um, 2)])     # mean U vector vs y
Ry = np.column_stack([avg(R, i) for i in range(6)])            # symmTensor xx,xy,xz,yy,yz,zz vs y
ywall = np.minimum(yc, 2*delta - yc)
Ly = np.minimum(0.41*ywall, 0.2*delta)                        # DFSEM length scale ~ min(kappa*y, 0.2 delta)
Ly = np.clip(Ly, 0.02*delta, None)

out = Path(args.out) if args.out else (case.parent/f"inflow_dfsem_{case.name}")
bd = out/"constant"/"boundaryData"/args.patch
(bd/"0").mkdir(parents=True, exist_ok=True)

def write_list(fn, rows, comps):
    with open(fn, "w") as f:
        f.write(f"{len(rows)}\n(\n")
        for r in rows:
            f.write("(" + " ".join(f"{v:.8g}" for v in (r if comps>1 else [r])) + ")\n" if comps>1 else f"{r:.8g}\n")
        f.write(")\n")
# points at (0, y, 0)
with open(bd/"points","w") as f:
    f.write(f"{len(yc)}\n(\n" + "\n".join(f"(0 {y:.8g} 0)" for y in yc) + "\n)\n")
write_list(bd/"0"/"U", Uy, 3)
write_list(bd/"0"/"R", Ry, 6)
with open(bd/"0"/"L","w") as f:
    f.write(f"{len(Ly)}\n(\n" + "\n".join(f"{v:.8g}" for v in Ly) + "\n)\n")

Ubar = float(Uy[:,0].mean())
bc = f"""// STEP 4: synthetic-turbulence inlet for {args.patch}, profiles from developed {case.name}/{t}
    {args.patch}
    {{
        type            turbulentDFSEMInlet;
        delta           {delta:.6g};
        nCellPerEddy    5;
        U               {{ type mappedFile; }}   // OF2512: U/R/L are PatchFunction1 mappedFile (verified by feature_smoke.py)
        R               {{ type mappedFile; }}
        L               {{ type mappedFile; }}
        value           uniform ({Ubar:.6g} 0 0);
    }}
// requires: constant/boundaryData/{args.patch}/points + 0/{{U,R,L}} (written here)
// NOTE: for a 2D geometry inlet, regenerate boundaryData as a PLANE with dfsem_remap_2d.py (1D line fails)."""
(out/"inlet_BC_snippet.txt").write_text(bc+"\n")

print(f"=== STEP 4 DFSEM inflow from {case.name}/{t} ===")
print(f"  delta={delta:.5g} m  Ubar={Ubar:.4g} m/s  profile points={len(yc)}")
print(f"  peak u'+ ~ {math.sqrt(max(Ry[:,0]))/math.sqrt(abs(0)+1e-9) if False else math.sqrt(max(Ry[:,0])):.4g} (raw <u'u'>^0.5 m/s)")
print(f"  wrote boundaryData -> {bd}")
print(f"  wrote inlet BC snippet -> {out/'inlet_BC_snippet.txt'}")
print(f"  -> drop boundaryData into the geometry tier's case; paste the BC into 0/U. One precursor feeds all tiers.")
