#!/usr/bin/env python3
"""Precursor library -- bank developed periodic-channel runs by Reynolds number so any geometry tier can pull a
MATCHED-Re turbulence inflow (the gap that made C11 silently run at the wrong Re). Pairs with inlet_contract.py.

  register <case> : scan a developed channel, compute (u_bar, delta, Re_bulk, Re_tau, Pr), append to index.json.
  match --re <Re> : return the banked precursor whose Re_bulk is closest to the tier target (+ % error).
  list            : show the bank.

index lives at precursor_library/index.json.
"""
import argparse, json, re, sys
from pathlib import Path
import numpy as np

WORK = Path(__file__).resolve().parent
LIB = WORK / "precursor_library"; LIB.mkdir(exist_ok=True)
INDEX = LIB / "index.json"

def rd(fn, kind):
    s = open(fn).read(); m = re.search(rf'internalField\s+nonuniform\s+List<{kind}>\s*\d+\s*\((.*?)\)\s*;', s, re.S)
    if not m: return None
    if kind == "scalar": return np.array([float(x) for x in m.group(1).split()])
    return np.array([[float(v) for v in n.split()] for n in re.findall(r'\(([^()]*)\)', m.group(1))])

def load(): return json.loads(INDEX.read_text()) if INDEX.exists() else []
def save(idx): INDEX.write_text(json.dumps(idx, indent=2) + "\n")

def scan(case, nu, pr, dpdx):
    case = Path(case)
    t = sorted([p.name for p in case.iterdir() if re.fullmatch(r"[0-9.]+", p.name)], key=float)[-1]
    Cy = rd(case/t/"Cy", "scalar"); Um = rd(case/t/"UMean", "vector")
    delta = float(Cy.max())/2; u_bar = float(Um[:, 0].mean())
    Re_bulk = u_bar*(2*delta)/nu
    if dpdx is None:  # try to read settled dpdx from the solver log
        for log in (case/"log.pimpleFoam",):
            if log.exists():
                ms = re.findall(r"pressure gradient = ([0-9.eE+-]+)", log.read_text())
                if ms: dpdx = float(ms[-1])
    u_tau = (dpdx*delta)**0.5 if dpdx else float("nan")
    Re_tau = u_tau*delta/nu if dpdx else float("nan")
    has_T = (case/t/"HMean").exists() or (case/t/"TMean").exists()
    return dict(name=case.name, path=str(case), time=t, u_bar=round(u_bar, 5), delta=delta,
                Re_bulk=round(Re_bulk, 1), Re_tau=round(Re_tau, 1) if dpdx else None,
                nu=nu, Pr=(pr if has_T else None), thermal=bool(has_T))

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("register"); r.add_argument("case"); r.add_argument("--nu", type=float, default=3.09278e-6)
    r.add_argument("--pr", type=float, default=5.0); r.add_argument("--dpdx", type=float, default=None)
    r.add_argument("--note", default="")
    m = sub.add_parser("match"); m.add_argument("--re", type=float, required=True); m.add_argument("--thermal", action="store_true")
    sub.add_parser("list")
    args = ap.parse_args()

    if args.cmd == "register":
        rec = scan(args.case, args.nu, args.pr, args.dpdx); rec["note"] = args.note
        idx = [e for e in load() if e["name"] != rec["name"]] + [rec]
        idx.sort(key=lambda e: e["Re_bulk"]); save(idx)
        print(f"registered {rec['name']}: Re_bulk={rec['Re_bulk']} Re_tau={rec['Re_tau']} "
              f"u_bar={rec['u_bar']} delta={rec['delta']} thermal={rec['thermal']}")
    elif args.cmd == "list":
        idx = load()
        print(f"=== precursor library ({len(idx)} banked) ===")
        for e in idx:
            print(f"  {e['name']:32s} Re_bulk={e['Re_bulk']:>8} Re_tau={str(e['Re_tau']):>6} "
                  f"u_bar={e['u_bar']:>6} thermal={e['thermal']} {e.get('note','')}")
    elif args.cmd == "match":
        idx = [e for e in load() if (e["thermal"] or not args.thermal)]
        if not idx: print("no precursors banked"); return 1
        best = min(idx, key=lambda e: abs(e["Re_bulk"]-args.re)/args.re)
        err = abs(best["Re_bulk"]-args.re)/args.re*100
        flag = "OK" if err <= 25 else "POOR (>25% -- run a new precursor at this Re)"
        print(f"target Re_bulk={args.re:.0f} -> best: {best['name']} (Re_bulk={best['Re_bulk']}, {err:.0f}% off) [{flag}]")
        print(f"  path={best['path']}  time={best['time']}  delta={best['delta']}  thermal={best['thermal']}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
