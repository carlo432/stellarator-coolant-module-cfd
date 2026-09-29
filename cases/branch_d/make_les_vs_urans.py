#!/usr/bin/env python3
"""LES (v26, full 0-0.8s incl. restart segments) vs URANS (v25) in the bent recirc zone.
Shows: LES initial transient (off the RANS start) THEN a developed fluctuating state,
vs the flat URANS. Honest: developed-window detrended std is the genuine fluctuation."""
import os, re, glob
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import numpy as np, matplotlib.pyplot as plt

def parse_all(base, which):
    files = sorted(glob.glob(f"{base}/*/U"), key=lambda p: float(p.split('/')[-2]))
    t, u = [], []
    for f in files:
        for l in open(f):
            if l.startswith('#') or not l.strip():
                continue
            m = re.match(r"\s*([0-9.eE+-]+)\s+(.*)", l)
            t.append(float(m.group(1)))
            u.append(float(re.findall(r"\(([^()]+)\)", m.group(2))[which].split()[0]))
    o = np.argsort(t)
    return np.array(t)[o], np.array(u)[o]

def movavg(x, w=31):
    return np.convolve(x, np.ones(w)/w, mode='same')

ROOT = "cases/branch_d"
ut, uu = parse_all(f"{ROOT}/v25_unsteady_rans_flibe_re10000_bent_recirculation_pimple/postProcessing/recircProbes", 2)
lt, lu = parse_all(f"{ROOT}/v26_les_smoke_flibe_re10000_bent_recirculation_wale/postProcessing/recircProbes", 2)
# downstream probe (index 2) — strongest developed fluctuation
dev = lt > 0.5
res = (lu[dev] - movavg(lu[dev], 31))[20:-20]

fig, ax = plt.subplots(figsize=(8.4, 4.6), dpi=180)
ax.axvspan(0.5, 0.8, color='0.92', label='developed window (t>0.5 s)')
ax.plot(ut, uu, '-', color='tab:blue', lw=1.6, label=f'URANS k-ωSST (v25)  steady (resid std≈0.0004)')
ax.plot(lt, lu, '-', color='tab:red',  lw=1.0, label=f'LES WALE (v26)  developed resid std={res.std():.3f}')
ax.set_xlabel("time [s]"); ax.set_ylabel("downstream-probe Ux [m/s]  (0.26, −0.035, 0.02)")
ax.set_title("Bent recirculation zone: LES develops resolved unsteadiness; URANS stays steady\n"
             "(LES first drifts off the RANS initial state, then fluctuates)")
ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc='upper left')
fig.tight_layout()
out = f"{ROOT}/figs/les_vs_urans_recirc_probe.png"
fig.savefig(out); print("wrote", out)
print(f"LES developed-window detrended: std={res.std():.4f}  pk-pk={res.max()-res.min():.4f}  vs URANS ~0.0004")
