#!/usr/bin/env python3
"""
C1 thermal certification verdict -- PRE-REGISTERED before the SGS-fix run landed.

Pre-registration (written 2026-06-25, BEFORE seeing the sgsthermal result), mirroring
how the O momentum cert was pre-registered:

  HYPOTHESIS: the +25% uniform H+ over-prediction in the no-SGS run is caused by the
  MISSING sub-grid thermal flux (scalarTransport used molecular D only). Adding it
  (alphaD=1/Pr, alphaDt=1/Prt_sgs) will move the LES H+ DOWN toward Kawamura DNS.

  PASS criteria (all must hold) for the SGS-fix run vs Kawamura Pr=5, Re_tau=180 DNS:
   P1 (PRIMARY)  mean |H+ error| <= 10% across the log layer, y+ in [30,150].
   P2 (sublayer) near-wall slope ratio T+/(Pr*y+) in [0.95,1.05] at y+<2 (T_tau sanity).
   P3 (center)   centerline H+ within +/-12% of DNS (~44.6; allow a touch low for Re_tau 167<180).
   P4 (Pr_t)     resolved Pr_t in [0.8,1.05] median over y+ in [20,100].
  DIRECTION   the SGS-fix H+ must be LOWER than the no-SGS H+ at the centerline
              (confirms the SGS flux adds mixing in the right direction).

  Outcomes:
   - all P1-P4 hold              -> CERTIFIED (thermal analog of the O cert).
   - direction holds, P1 misses  -> "SGS flux is necessary but not sufficient":
                                    points to residual near-wall thermal under-resolution
                                    and/or Prt_sgs tuning -> quantifies the HPC-1 requirement.
   - LES now BELOW DNS           -> Prt_sgs=0.9 over-mixes -> raise Prt_sgs.
"""
import argparse, csv, json, math
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORK = Path(__file__).resolve().parent
DIG = WORK / "kawamura_pr5_dns_digitized"
PR = 5.0


def col(path, xk, yk):
    xs, ys = [], []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try: x, y = float(r[xk]), float(r[yk])
            except (ValueError, KeyError): continue
            if math.isfinite(x) and math.isfinite(y): xs.append(x); ys.append(y)
    return xs, ys


def interp(xs, ys, xq):
    if xq <= xs[0]: return ys[0]
    if xq >= xs[-1]: return ys[-1]
    for i in range(1, len(xs)):
        if xs[i-1] <= xq <= xs[i]:
            t = (xq-xs[i-1])/(xs[i]-xs[i-1]); return ys[i-1]+t*(ys[i]-ys[i-1])
    return ys[-1]


def median(v): v=sorted(v); n=len(v); return float("nan") if not v else (v[n//2] if n%2 else 0.5*(v[n//2-1]+v[n//2]))

ap = argparse.ArgumentParser()
ap.add_argument("--nosgs", required=True)
ap.add_argument("--sgs", required=True)
ap.add_argument("--json")
ap.add_argument("--plot")
args = ap.parse_args()

dy, dT = col(DIG/"kawamura_fig2_mean_temperature.csv", "y_plus", "Theta_plus")
dyp, dprt = col(DIG/"kawamura_fig9_turbulent_prandtl.csv", "y_plus", "Pr_t")
ny, nT = col(args.nosgs, "y_plus", "T_plus")
sy, sT = col(args.sgs, "y_plus", "T_plus")
syp, sprt = col(args.sgs, "y_plus", "Pr_t_resolved")

# --- pre-registered metrics on the SGS-fix run ---
loglayer = [y for y in (30, 50, 75, 100, 150)]
errs = [abs(interp(sy, sT, y) - interp(dy, dT, y)) / interp(dy, dT, y) for y in loglayer]
P1_val = max(errs) * 100
P1 = P1_val <= 10.0
slope_ratio = interp(sy, sT, 1.79) / (PR * 1.79)
P2 = 0.95 <= slope_ratio <= 1.05
center_les, center_dns = sT[-1], dT[-1]
P3_val = (center_les - center_dns) / center_dns * 100
P3 = abs(P3_val) <= 12.0
prt_band = [p for y, p in zip(syp, sprt) if 20 <= y <= 100 and math.isfinite(p)]
P4_val = median(prt_band)
P4 = 0.8 <= P4_val <= 1.05
direction = center_les < nT[-1]

print("==================  C1 THERMAL CERT VERDICT (pre-registered)  ==================")
print(f"  no-SGS centerline H+ = {nT[-1]:.1f}   SGS-fix centerline H+ = {center_les:.1f}   DNS ~ {center_dns:.1f}")
print(f"  DIRECTION (SGS lowers H+ vs no-SGS): {'YES' if direction else 'NO'}")
print(f"  P1 log-layer max|err| = {P1_val:5.1f}%  (<=10%)   -> {'PASS' if P1 else 'FAIL'}")
print(f"  P2 sublayer slope T+/(Pr y+)@1.79 = {slope_ratio:.3f}  [0.95,1.05] -> {'PASS' if P2 else 'FAIL'}")
print(f"  P3 centerline err = {P3_val:+5.1f}%  (+/-12%)      -> {'PASS' if P3 else 'FAIL'}")
print(f"  P4 Pr_t median(y+20-100) = {P4_val:.3f}  [0.8,1.05]  -> {'PASS' if P4 else 'FAIL'}")
verdict = "CERTIFIED" if (P1 and P2 and P3 and P4) else ("SGS-NECESSARY-NOT-SUFFICIENT" if direction else "OVER/UNDER-MIXED")
print(f"  ------> VERDICT: {verdict}")
print("================================================================================")

result = {
    "nosgs_profile": str(Path(args.nosgs)),
    "sgs_profile": str(Path(args.sgs)),
    "criteria": {
        "P1": {"value_percent": P1_val, "limit_percent": 10.0, "passed": P1},
        "P2": {"value": slope_ratio, "limits": [0.95, 1.05], "passed": P2},
        "P3": {"value_percent": P3_val, "limits_percent": [-12.0, 12.0], "passed": P3},
        "P4": {"value": P4_val, "limits": [0.8, 1.05], "passed": P4},
    },
    "direction_passed": direction,
    "centerline": {"nosgs": nT[-1], "sgs": center_les, "dns": center_dns},
    "verdict": verdict,
}
if args.json:
    json_path = Path(args.json)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print("wrote", json_path)

# 3-way overlay
import numpy as np
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
a1.semilogx(dy, dT, "o", mfc="none", color="k", label="Kawamura DNS Pr=5")
a1.semilogx(ny, nT, "-^", ms=3, color="tab:red", label="LES no-SGS (molecular D)")
a1.semilogx(sy, sT, "-s", ms=3, color="tab:green", label="LES +SGS thermal flux")
yk = np.logspace(0, math.log10(6), 40); a1.semilogx(yk, [PR*y for y in yk], "--", color="tab:blue", lw=1, label=r"$\Theta^+=Pr\,y^+$")
a1.set_xlabel(r"$y^+$"); a1.set_ylabel(r"$H^+$"); a1.set_title(f"C1 mean temp -- {verdict}")
a1.set_ylim(0, max(nT)*1.05); a1.grid(True, which="both", alpha=0.3); a1.legend(fontsize=8)
a2.semilogx(dyp, dprt, "o", mfc="none", color="k", label="DNS Pr_t")
a2.semilogx(syp, sprt, "-s", ms=3, color="tab:green", label="LES +SGS Pr_t")
a2.axhline(0.85, ls=":", color="gray"); a2.set_ylim(0, 1.6); a2.set_xlabel(r"$y^+$"); a2.set_ylabel(r"$Pr_t$")
a2.set_title("Turbulent Prandtl"); a2.grid(True, which="both", alpha=0.3); a2.legend(fontsize=8)
out = Path(args.plot) if args.plot else DIG/"c1_thermal_verdict_3way.png"
out.parent.mkdir(parents=True, exist_ok=True)
fig.tight_layout(); fig.savefig(out, dpi=130)
print("wrote", out)
