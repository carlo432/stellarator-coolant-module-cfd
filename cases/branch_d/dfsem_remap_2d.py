#!/usr/bin/env python3
"""STEP 4 fix-up: remap a 1D periodic-channel precursor profile onto a 2D geometry-tier inlet PLANE.

hpc_inflow_dfsem.py emits a 1D profile (points on a line (0,y,0)) -- correct for a channel, but a geometry
tier's inlet is a 2D patch, and turbulentDFSEMInlet needs the boundaryData points to span a PLANE (it builds a
plane normal). This reads the meshed case's actual inlet face centres, maps the channel wall-normal profile onto
the inlet's wall-normal axis (the SHORTER in-plane extent = the no-slip gap), tiles it across the span, and writes
boundaryData at the exact face centres so mapMethod nearestCell is 1:1.

Usage: dfsem_remap_2d.py <geometry_case> <precursor_channel_case> [--patch inlet] [--time latestTime]
"""
import argparse, re, sys
from pathlib import Path
import numpy as np

def read_internal(fn, kind):
    s = open(fn).read()
    m = re.search(rf'internalField\s+nonuniform\s+List<{kind}>\s*\d+\s*\((.*?)\)\s*;', s, re.S)
    if kind == "scalar":
        return np.array([float(x) for x in m.group(1).split()])
    return np.array([[float(v) for v in n.split()] for n in re.findall(r'\(([^()]*)\)', m.group(1))])

def channel_profile(case, t):
    Cy = read_internal(case/t/"Cy", "scalar")
    Um = read_internal(case/t/"UMean", "vector")
    R  = read_internal(case/t/"UPrime2Mean", "symmTensor")
    ylv = np.unique(np.round(Cy, 9))
    sel = lambda a, c: np.array([a[np.isclose(Cy, y), c].mean() for y in ylv])
    yc = np.array([Cy[np.isclose(Cy, y)].mean() for y in ylv])
    U  = np.column_stack([sel(Um, i) for i in range(3)])
    Rr = np.column_stack([sel(R, i) for i in range(6)])
    delta = yc.max()/2
    ywall = np.minimum(yc, 2*delta - yc)
    L  = np.clip(np.minimum(0.41*ywall, 0.2*delta), 0.02*delta, None)
    return yc, U, Rr, L, delta

def inlet_centres(case, patch):
    bf = (case/"constant/polyMesh/boundary").read_text()
    m = re.search(rf"{patch}\s*\{{[^}}]*?nFaces\s+(\d+);[^}}]*?startFace\s+(\d+);", bf, re.S)
    nF, sF = int(m.group(1)), int(m.group(2))
    pts = (case/"constant/polyMesh/points").read_text()
    parr = np.array([[float(v) for v in t.split()] for t in re.findall(r"\(([-0-9.eE ]+)\)", pts.split("(",1)[1])])
    fc = (case/"constant/polyMesh/faces").read_text()
    faces = re.findall(r"\d+\(([\d ]+)\)", fc)
    cent = np.array([parr[[int(i) for i in f.split()]].mean(0) for f in faces[sF:sF+nF]])
    return cent

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case"); ap.add_argument("precursor")
    ap.add_argument("--patch", default="inlet"); ap.add_argument("--time", default=None)
    args = ap.parse_args()
    case = Path(args.case); pre = Path(args.precursor)
    t = args.time or sorted([p.name for p in pre.iterdir() if re.fullmatch(r"[0-9.]+", p.name)], key=float)[-1]

    yc, U, R, L, delta = channel_profile(pre, t)
    cent = inlet_centres(case, args.patch)
    # pick wall-normal axis = the in-plane axis with the SHORTER extent (the no-slip gap); span = the longer
    spans = [np.ptp(cent[:, a]) for a in range(3)]
    flow_ax = int(np.argmin([len(np.unique(np.round(cent[:, a], 7))) for a in range(3)]))  # constant axis = flow normal
    inplane = [a for a in range(3) if a != flow_ax]
    wn = inplane[0] if spans[inplane[0]] <= spans[inplane[1]] else inplane[1]   # wall-normal = shorter
    w = cent[:, wn]
    # map inlet wall-normal coord linearly to channel yc in [0, 2 delta]
    chan = (w - w.min())/(w.max()-w.min()) * (2*delta)
    def interp_cols(arr):
        return np.column_stack([np.interp(chan, yc, arr[:, c]) for c in range(arr.shape[1])])
    Uf = interp_cols(U); Rf = interp_cols(R); Lf = np.interp(chan, yc, L)

    bd = case/"constant"/"boundaryData"/args.patch
    (bd/"0").mkdir(parents=True, exist_ok=True)
    with open(bd/"points","w") as f:
        f.write(f"{len(cent)}\n(\n" + "\n".join(f"({p[0]:.8g} {p[1]:.8g} {p[2]:.8g})" for p in cent) + "\n)\n")
    def wr(fn, arr, nc):
        with open(fn,"w") as f:
            f.write(f"{len(arr)}\n(\n")
            for r in (arr if nc>1 else arr.reshape(-1,1)):
                f.write(("("+" ".join(f"{v:.8g}" for v in r)+")\n") if nc>1 else f"{r[0]:.8g}\n")
            f.write(")\n")
    wr(bd/"0"/"U", Uf, 3); wr(bd/"0"/"R", Rf, 6); wr(bd/"0"/"L", Lf, 1)
    print(f"=== DFSEM 2D remap: {pre.name}/{t} -> {case.name}/{args.patch} ===")
    print(f"  inlet faces={len(cent)}  flow_axis={'xyz'[flow_ax]}  wall_normal={'xyz'[wn]} (gap={spans[wn]:.4g})  delta_chan={delta:.4g}")
    print(f"  Ubar_inlet={Uf[:,0].mean():.4g}  peak<u'u'>^.5={np.sqrt(Rf[:,0].max()):.4g}")
    print(f"  wrote 2D plane boundaryData -> {bd}")

if __name__ == "__main__":
    main()
