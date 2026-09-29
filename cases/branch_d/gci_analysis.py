#!/usr/bin/env python3
"""3-grid Richardson/GCI (ASME) on peak heated-wall dT for the outlet-offset case.
Usage: python3 gci_analysis.py <peak_fine>   (coarse/medium baked in; fine from arg)"""
import sys, math
# (cells, peak dT) ; 1=fine, 2=medium, 3=coarse
N3, p3 = 21248, 94.43     # coarse
N2, p2 = 73128, 77.00     # medium (v23)
N1 = 218080
p1 = float(sys.argv[1]) if len(sys.argv) > 1 else None

def gci(N1, N2, N3, f1, f2, f3, Fs=1.25):
    r21 = (N1 / N2) ** (1/3); r32 = (N2 / N3) ** (1/3)
    e32 = f2 - f3; e21 = f1 - f2
    s = math.copysign(1, e32 / e21) if e21 != 0 else 1
    p = 1.5
    for _ in range(100):
        q = math.log((r21**p - s) / (r32**p - s))
        p = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
    f_ext = (r21**p * f1 - f2) / (r21**p - 1)
    ea21 = abs((f1 - f2) / f1)
    gci21 = Fs * ea21 / (r21**p - 1)
    return dict(r21=r21, r32=r32, p=p, f_ext=f_ext, gci21_pct=gci21*100,
                err_ext_pct=abs((f_ext - f1) / f_ext)*100)

print(f"coarse {N3:>7d} cells: peak dT = {p3:.2f} K")
print(f"medium {N2:>7d} cells: peak dT = {p2:.2f} K")
if p1 is None:
    print("fine: <pending>  (rerun with: python3 gci_analysis.py <peak_fine>)")
else:
    print(f"fine   {N1:>7d} cells: peak dT = {p1:.2f} K")
    g = gci(N1, N2, N3, p1, p2, p3)
    print(f"\nrefinement ratios r21={g['r21']:.3f} r32={g['r32']:.3f}")
    print(f"observed order p = {g['p']:.2f}")
    print(f"Richardson-extrapolated peak dT (h->0) = {g['f_ext']:.2f} K")
    print(f"GCI_fine = {g['gci21_pct']:.2f} %  (uncertainty band on the fine peak dT)")
    print(f"fine-vs-extrapolated error = {g['err_ext_pct']:.2f} %")
    print(f"\n=> medium(73k)={p2:.1f}, fine={p1:.1f}, extrapolated={g['f_ext']:.1f} K. "
          f"Ordering check: bend 59.9 < straight 63.5 < offset/BFS ~77+/-GCI -- holds if band < ~13 K.")
