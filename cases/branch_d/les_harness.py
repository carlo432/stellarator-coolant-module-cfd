#!/usr/bin/env python3
"""Spine M: LES verification harness. Reads a converged time dir with fieldAverage output (UMean/UPrime2Mean)
+ cell centres (Cy), plane-averages into a wall-unit mean-velocity profile, computes u_tau/Re_tau, and checks
it against the law of the wall (viscous sublayer U+=y+ and log law U+=(1/kappa)ln(y+)+B). Reusable on any
channel LES (L now, HPC-1 later)."""
import re, sys, math, statistics
from pathlib import Path
import numpy as np

gate_only = "--scalar-gate-only" in sys.argv
positional = [arg for arg in sys.argv[1:] if not arg.startswith("--")]
case = positional[0] if positional else "les_channel_wm"
t    = positional[1] if len(positional)>1 else "2.09197"
nu   = 3.09278e-6
delta= 0.01
dpdx = float(positional[2]) if len(positional)>2 else 0.35392258   # kinematic driving pressure gradient [m/s^2] (per-case; from solver log)
kappa, B = 0.41, 5.2

case_path = Path(case)
control_dict = case_path / "system" / "controlDict"
control_text = control_dict.read_text() if control_dict.exists() else ""
contract_path = case_path / "constant" / "sourceDrivenScalarGate.json"
has_source = "scalarSemiImplicitSource" in control_text
if has_source or contract_path.exists() or gate_only:
    if not contract_path.exists():
        print(
            "STATISTICS REFUSED: source-driven scalar has no "
            "constant/sourceDrivenScalarGate.json contract."
        )
        raise SystemExit(2)
    from scalar_convergence_gate import evaluate_contract, print_report

    gate_result = evaluate_contract(case_path, contract_path)
    print_report(gate_result)
    if not gate_result["statistics_allowed"]:
        raise SystemExit(2)
    if gate_only:
        raise SystemExit(0)

def read_scalar(fn):
    s=open(fn).read()
    m=re.search(r'internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)\s*;', s, re.S)
    return np.array([float(x) for x in m.group(1).split()])

def read_vector(fn):
    s=open(fn).read()
    m=re.search(r'internalField\s+nonuniform\s+List<vector>\s*\d+\s*\((.*?)\)\s*;', s, re.S)
    nums=re.findall(r'\(([^()]*)\)', m.group(1))
    return np.array([[float(v) for v in n.split()] for n in nums])

def read_symmtensor(fn):   # 6 comp per entry: (xx xy xz yy yz zz)
    s=open(fn).read()
    m=re.search(r'internalField\s+nonuniform\s+List<symmTensor>\s*\d+\s*\((.*?)\)\s*;', s, re.S)
    nums=re.findall(r'\(([^()]*)\)', m.group(1))
    return np.array([[float(v) for v in n.split()] for n in nums])

Cy   = read_scalar(f"{case}/{t}/Cy")
Um   = read_vector(f"{case}/{t}/UMean")          # mean velocity vectors
Umx  = Um[:,0]

# group cells by y-level (64 levels, ~1024 cells each)
ylevels = np.unique(np.round(Cy, 9))
yc   = np.array([Cy[np.isclose(Cy,y)].mean() for y in ylevels])
Uprof= np.array([Umx[np.isclose(Cy,y)].mean() for y in ylevels])

# u_tau two ways: (1) pressure-gradient balance, (2) near-wall mean gradient
utau_pg = math.sqrt(abs(dpdx)*delta)
# wall gradient from the lowest cell (y measured from bottom wall)
ybot = yc.min(); Ubot = Uprof[np.argmin(yc)]
utau_wg = math.sqrt(nu*Ubot/ybot)
Re_tau  = utau_pg*delta/nu

# fold both walls: wall distance = min(y, 2delta - y)
ywall = np.minimum(yc, 2*delta - yc)
yp    = ywall*utau_pg/nu
Up    = Uprof/utau_pg
# keep lower half (unique wall-distance), sort
half  = yc <= delta
order = np.argsort(yp[half])
ypH, UpH = yp[half][order], Up[half][order]

print(f"=== LES harness: {case}/{t} ===")
print(f"u_tau (pressure-grad) = {utau_pg:.5f} m/s   u_tau (wall-grad) = {utau_wg:.5f} m/s   ratio={utau_wg/utau_pg:.2f}")
print(f"Re_tau = {Re_tau:.1f}   (nearest MKM DNS case: {'180' if Re_tau<290 else '395'})")
print(f"first-cell y+ = {ypH[0]:.2f}   (wall-resolved if <~1-2)")
# log-law check in the log region (30 < y+ < ~0.3 Re_tau)
mask=(ypH>30)&(ypH<0.4*Re_tau)
if mask.sum()>=2:
    loglaw=(1/kappa)*np.log(ypH[mask])+B
    err=np.mean(np.abs(UpH[mask]-loglaw)/loglaw)*100
    print(f"log-law region pts={mask.sum()}  mean |U+ - loglaw|/loglaw = {err:.1f}%")
else:
    print(f"log-law region: too few pts (Re_tau low); compare full profile to MKM-180 instead")

# ---- figure ----
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(7.2,5.4))
ax.semilogx(ypH, UpH, 'o-', ms=4, color='#C44E52', label=f'wall-modeled LES (Re_τ≈{Re_tau:.0f})')
ysub=np.logspace(-1,1.1,50); ax.semilogx(ysub, ysub, 'k--', lw=1, label='viscous sublayer U⁺=y⁺')
ylog=np.logspace(math.log10(8),math.log10(max(ypH)),50)
ax.semilogx(ylog, (1/kappa)*np.log(ylog)+B, 'b-.', lw=1.2, label=f'log law (1/{kappa})ln y⁺+{B}')
ax.set_xlabel('y⁺'); ax.set_ylabel('U⁺')
ax.set_title(f'Spine M: LES mean-velocity profile vs law of the wall\n{case} t={t}, ~20 flow-through average, FLiBe channel')
ax.legend(loc='upper left', fontsize=9); ax.grid(alpha=0.3, which='both')
ax.set_xlim(0.3, max(ypH)*1.3)
plt.tight_layout(); plt.savefig(f"{case}/../figs/m_les_wall_profile.png", dpi=150)
print(f"wrote figs/m_les_wall_profile.png")

# ============ resolved Reynolds-stress / turbulence-intensity profiles ============
# UPrime2Mean (symmTensor) = <u'_i u'_j> resolved. Normal stresses -> RMS intensities in wall units.
R = read_symmtensor(f"{case}/{t}/UPrime2Mean")     # cols: xx,xy,xz,yy,yz,zz
def yprof(col): return np.array([R[np.isclose(Cy,y),col].mean() for y in ylevels])
Ruu,Rvv,Rww = yprof(0),yprof(3),yprof(5)
urms=np.sqrt(np.clip(Ruu,0,None))/utau_pg          # u'+  (streamwise)
vrms=np.sqrt(np.clip(Rvv,0,None))/utau_pg          # v'+  (wall-normal)
wrms=np.sqrt(np.clip(Rww,0,None))/utau_pg          # w'+  (spanwise)
tke =0.5*(Ruu+Rvv+Rww)/utau_pg**2                  # resolved TKE+
# lower half, sorted by y+
uH,vH,wH,kH = urms[half][order],vrms[half][order],wrms[half][order],tke[half][order]
ipk=np.argmax(uH)
print(f"resolved u'+ peak = {uH[ipk]:.2f} at y+ = {ypH[ipk]:.1f}   (DNS channel ~2.6-2.8 at y+~12-15; lower=under-resolved)")
print(f"resolved TKE+ peak = {kH.max():.2f}")

fig2,ax2=plt.subplots(figsize=(7.2,5.4))
ax2.semilogx(ypH,uH,'o-',ms=3,color='#C44E52',label="u'⁺ (streamwise)")
ax2.semilogx(ypH,vH,'s-',ms=3,color='#4C72B0',label="v'⁺ (wall-normal)")
ax2.semilogx(ypH,wH,'^-',ms=3,color='#55A868',label="w'⁺ (spanwise)")
ax2.axhspan(2.6,2.8,color='grey',alpha=0.18); ax2.text(1.1,2.85,'DNS u\'⁺ peak band ~2.6-2.8',fontsize=8,color='#555')
ax2.set_xlabel('y⁺'); ax2.set_ylabel("resolved RMS velocity  (wall units)")
ax2.set_title(f'Spine M: resolved turbulence intensities vs y⁺\n{case} t={t}, Re_τ≈{Re_tau:.0f}')
ax2.legend(loc='upper right',fontsize=9); ax2.grid(alpha=0.3,which='both'); ax2.set_xlim(0.3,max(ypH)*1.3)
plt.tight_layout(); plt.savefig(f"{case}/../figs/m_les_reystress_{case.split('_')[-1]}.png", dpi=150)
print(f"wrote figs/m_les_reystress_{case.split('_')[-1]}.png")
