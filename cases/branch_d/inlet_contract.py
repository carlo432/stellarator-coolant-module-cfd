#!/usr/bin/env python3
"""BUCKET B -- physics/interface de-risk. BEFORE wiring a periodic precursor into a geometry tier's inlet,
assert the three things that silently bit us on C11:

  (a) Re MATCH      : precursor bulk Reynolds number == tier target Re (DFSEM imposes the PRECURSOR's mean;
                      C8 ran at u_bar~2.0 while the tier wanted 0.86 -> the tier silently ran at the wrong Re).
  (b) DIMENSIONALITY : the inlet boundaryData must span a 2D PLANE, not a 1D line (turbulentDFSEMInlet builds a
                      plane normal; a channel's 1D profile is a line and fails: 'points on a single line').
  (c) delta <-> gap  : the precursor half-height delta must map to the inlet wall-normal gap (so the profile
                      lands on the right length scale).

Read-only. Exit non-zero if any HARD check (a/b) fails.

Usage: inlet_contract.py <geometry_case> <precursor_channel_case> [--patch inlet] [--u-tier <m/s>]
                         [--nu 3.09278e-6] [--re-tol 0.25]

If mapped boundaryData exists, the contract measures the delivered bulk speed
from `constant/boundaryData/<patch>/0/U` and uses that for the Re check.  The
`--u-tier` argument is only a fallback for cases without mapped inlet data.
"""
import argparse, re, sys
from pathlib import Path
import numpy as np

def read_internal(fn, kind):
    s = open(fn).read()
    m = re.search(rf'internalField\s+nonuniform\s+List<{kind}>\s*\d+\s*\((.*?)\)\s*;', s, re.S)
    if not m: return None
    if kind == "scalar":
        return np.array([float(x) for x in m.group(1).split()])
    return np.array([[float(v) for v in n.split()] for n in re.findall(r'\(([^()]*)\)', m.group(1))])

def precursor_state(pre, nu):
    t = sorted([p.name for p in pre.iterdir() if re.fullmatch(r"[0-9.]+", p.name)], key=float)[-1]
    Cy = read_internal(pre/t/"Cy", "scalar"); Um = read_internal(pre/t/"UMean", "vector")
    delta = Cy.max()/2
    Ub = float(Um[:, 0].mean())             # bulk ~ volume-mean streamwise (uniform channel)
    Re = Ub*(2*delta)/nu
    return t, Ub, delta, Re

def inlet_centres(case, patch):
    bf = (case/"constant/polyMesh/boundary").read_text()
    m = re.search(rf"{patch}\s*\{{[^}}]*?nFaces\s+(\d+);[^}}]*?startFace\s+(\d+);", bf, re.S)
    nF, sF = int(m.group(1)), int(m.group(2))
    pts = (case/"constant/polyMesh/points").read_text()
    parr = np.array([[float(v) for v in t.split()] for t in re.findall(r"\(([-0-9.eE ]+)\)", pts.split("(",1)[1])])
    faces = re.findall(r"\d+\(([\d ]+)\)", (case/"constant/polyMesh/faces").read_text())
    return np.array([parr[[int(i) for i in f.split()]].mean(0) for f in faces[sF:sF+nF]])

def read_vector_list(path):
    if not path.exists():
        return None
    text = path.read_text()
    return np.array([[float(v) for v in n.split()] for n in re.findall(r'\(([^()]*)\)', text)])

def boundary_data_bulk(case, patch):
    vecs = read_vector_list(case/"constant"/"boundaryData"/patch/"0"/"U")
    if vecs is None or len(vecs) == 0:
        return None
    return float(np.linalg.norm(vecs, axis=1).mean())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case"); ap.add_argument("precursor")
    ap.add_argument("--patch", default="inlet"); ap.add_argument("--nu", type=float, default=3.09278e-6)
    ap.add_argument("--u-tier", type=float, default=None, help="tier target bulk U [m/s]; else read inlet/internalField")
    ap.add_argument("--re-tol", type=float, default=0.25)
    args = ap.parse_args()
    case = Path(args.case); pre = Path(args.precursor)

    t, Ub, delta, Re_pre = precursor_state(pre, args.nu)
    cent = inlet_centres(case, args.patch)
    spans = [np.ptp(cent[:, a]) for a in range(3)]
    nuniq = [len(np.unique(np.round(cent[:, a], 7))) for a in range(3)]
    flow_ax = int(np.argmin(nuniq))
    inplane = [a for a in range(3) if a != flow_ax]
    gap = min(spans[inplane[0]], spans[inplane[1]])         # wall-normal gap = shorter in-plane extent

    # tier target U: mapped boundaryData delivered bulk wins; arg/internalField are fallbacks.
    u_source = "boundaryData mean|U|"
    u_tier = boundary_data_bulk(case, args.patch)
    if u_tier is None:
        u_source = "--u-tier"
        u_tier = args.u_tier
    if u_tier is None:
        u_source = "0/U internalField"
        m = re.search(r"internalField\s+uniform\s+\(([^)]+)\)", (case/"0/U").read_text())
        u_tier = float(m.group(1).split()[0]) if m else Ub
    if u_tier is None:
        u_source = "precursor fallback"
        u_tier = Ub
    Re_tier = u_tier*gap/args.nu

    # boundaryData dimensionality (if present)
    bdp = case/"constant"/"boundaryData"/args.patch/"points"
    dim = None
    if bdp.exists():
        pts = np.array([[float(v) for v in s.split()] for s in re.findall(r"\(([-0-9.eE ]+)\)", bdp.read_text())])
        varia = sum(1 for a in range(3) if len(np.unique(np.round(pts[:, a], 7))) > 1)
        dim = "plane(2D)" if varia >= 2 else "line(1D)"

    checks = []
    re_err = abs(Re_pre - Re_tier)/Re_tier
    checks.append(("a Re match", re_err <= args.re_tol, True,
                   f"Re_pre={Re_pre:.0f} (Ub={Ub:.3g}) vs Re_tier={Re_tier:.0f} (U={u_tier:.3g}, source={u_source}) -> {re_err*100:.0f}% off"))
    if dim is None:
        checks.append(("b dimensionality", True, False, "no boundaryData yet (run dfsem_remap_2d.py to build the PLANE)"))
    else:
        checks.append(("b dimensionality", dim == "plane(2D)", True, f"inlet boundaryData is a {dim} ({len(cent)} faces)"))
    ratio = (2*delta)/gap
    checks.append(("c delta<->gap", 0.5 <= ratio <= 2.0, False, f"2*delta_pre={2*delta:.4g} vs inlet gap={gap:.4g} -> ratio {ratio:.2f}"))

    print(f"=== inlet contract: {pre.name}/{t}  ->  {case.name}/{args.patch} ===")
    hardfail = 0
    for name, ok, hard, detail in checks:
        mark = "[ OK ]" if ok else ("[FAIL]" if hard else "[WARN]")
        print(f"  {mark} {name}: {detail}")
        if hard and not ok: hardfail += 1
    print(f"--- {hardfail} hard failure(s) ---" if hardfail else "--- contract satisfied (hard checks) ---")
    return 1 if hardfail else 0

if __name__ == "__main__":
    sys.exit(main())
