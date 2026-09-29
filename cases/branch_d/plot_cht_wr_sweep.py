#!/usr/bin/env python3
"""Plot the wall-resolved CHT velocity sweep: the trustworthy-metric version of the
flushing-speed law. Panel A: interface wall dT vs inlet U (fit + U^-0.8 reference).
Panel B: Nu_CHT vs Re vs Dittus-Boelter & Gnielinski correlations.
Usage: python3 plot_cht_wr_sweep.py [results.txt]"""
import sys, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

f = sys.argv[1] if len(sys.argv) > 1 else "cht_wr_sweep_results.txt"
rows = [l.split() for l in open(f) if l.strip() and not l.startswith("#")]
U  = [float(r[0]) for r in rows]
Re = [float(r[1]) for r in rows]
dT = [float(r[3]) for r in rows]
Nu = [float(r[4]) for r in rows]

# power-law fit dT = a * U^n
lx = [math.log(u) for u in U]; ly = [math.log(d) for d in dT]
n  = len(U); sx=sum(lx); sy=sum(ly); sxx=sum(x*x for x in lx); sxy=sum(a*b for a,b in zip(lx,ly))
slope = (n*sxy - sx*sy)/(n*sxx - sx*sx)
inter = (sy - slope*sx)/n
a = math.exp(inter)

Pr, Dh, kf, q = 14.4, 0.02667, 1.0, 100000.0
def nu_db(re): return 0.023*re**0.8*Pr**0.4
def nu_g(re):
    fr=(0.790*math.log(re)-1.64)**-2
    return (fr/8)*(re-1000)*Pr/(1+12.7*math.sqrt(fr/8)*(Pr**(2/3)-1))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Panel A
ax1.plot(U, dT, "o", ms=9, color="#c0392b", label="wall-resolved CHT", zorder=5)
uu = [min(U)*0.95 + i*(max(U)*1.05-min(U)*0.95)/100 for i in range(101)]
ax1.plot(uu, [a*u**slope for u in uu], "-", color="#c0392b",
         label=f"fit  ΔT ∝ U^{slope:.2f}")
a08 = dT[2]/U[2]**-0.8
ax1.plot(uu, [a08*u**-0.8 for u in uu], "--", color="#2c3e50",
         label="Dittus-Boelter  ΔT ∝ U^-0.8")
ax1.set_xlabel("inlet velocity U (m/s)"); ax1.set_ylabel("interface wall ΔT (K)")
ax1.set_title("Wall superheat vs flushing speed\n(trustworthy CHT metric, 100 kW/m²)")
ax1.legend(); ax1.grid(alpha=0.3)

# Panel B
ax2.plot(Re, Nu, "o", ms=9, color="#c0392b", label="Nu_CHT", zorder=5)
rr = [min(Re)*0.95 + i*(max(Re)*1.05-min(Re)*0.95)/100 for i in range(101)]
ax2.plot(rr, [nu_db(r) for r in rr], "--", color="#2c3e50", label="Dittus-Boelter")
ax2.plot(rr, [nu_g(r)  for r in rr], ":",  color="#16a085", lw=2, label="Gnielinski")
ax2.set_xlabel("Re"); ax2.set_ylabel("Nu"); ax2.set_title("CHT Nusselt vs correlations\n(Pr=14.4 FLiBe)")
ax2.legend(); ax2.grid(alpha=0.3)

plt.tight_layout()
out = "figs/cht_wr_velocity_sweep.png"
plt.savefig(out, dpi=140)
print(f"wrote {out}")
print(f"fitted exponent n = {slope:.3f} (target -0.8)")
print("mean |Nu_CHT - Gnielinski|/Gnielinski = "
      f"{sum(abs(Nu[i]-nu_g(Re[i]))/nu_g(Re[i]) for i in range(len(Re)))/len(Re)*100:.1f}%")
