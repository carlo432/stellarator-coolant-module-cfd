#!/usr/bin/env python3
"""FLUSHING-PRESCRIPTION CFD TEST: a locally contracted throat at the flux peak.

The flushing remedy (analyze_flushing_remedy.py) claims boiling is avoided by a 1.33x LOCAL
near-wall speed-up over a ~34 mm strip at the flux peak, sized from the DEVELOPED-duct law
dT_wall ~ U^-0.70. That law was measured on ducts that run at one speed everywhere. This case
DELIBERATELY STRESSES the local-similarity assumption behind it: a real contraction accelerates
the flow but also re-develops the boundary layer and can separate on the downstream expansion --
neither of which a developed-duct correlation knows about.

CONSTRUCTION (clean single-variable test): clone the straight plasma baseline
(geom_straight_pipeline_plasma) at its latest time and DEFORM ONLY THE COLD-SIDE (z<0) mesh
points to pinch the channel at the flux peak. Everything the wall sees is held fixed:
  - heated wall stays the z=+0.015 plane; its AREA is unchanged (deforming z<0 does not touch it);
  - the exact-flux plasma gradientExpr q''(s), s=pos().x(), is unchanged;
  - the OpenMC depth zones grow from the heated wall inward and their heated-side portion is
    unchanged; only the deep cold-side cells shrink (~3% total volumetric power at the throat,
    2nd-order for the wall hotspot -- noted, not corrected).
So the ONLY thing that changes is the near-wall flow: continuity forces U up by the target
factor at the throat. Baseline for comparison = geom_straight_pipeline_plasma itself.

Contraction: channel height 30 mm (z in [-15,+15] mm). Heated side z>0 fixed (thickness 15 mm);
cold side z<0 scaled by f(x) so throat height = 15*(1+f_center) mm. Speedup m = 2/(1+f_center)
=> f_center = 2/m - 1 (m=1.33 -> f_center=0.5 -> throat height 22.5 mm -> 1.333x). f(x) is a
smooth Gaussian bump of FWHM = strip width centered on the flux peak, so the shoulders are gentle
(low non-orthogonality) -- checkMesh is run and reported before any solve.

Usage: make_throat_pipeline_case.py [--speedup 1.33] [--strip-mm 34] [--end-time 4.0]
"""
from __future__ import annotations
import argparse, json, math, re, shutil, stat, subprocess
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--src", default="geom_straight_pipeline_plasma")
ap.add_argument("--case", default="geom_throat_pipeline_plasma")
ap.add_argument("--qmap", default="figs/plasma_qmap.json")
ap.add_argument("--speedup", type=float, default=1.33, help="target near-wall bulk speed-up at the throat")
ap.add_argument("--strip-mm", type=float, default=34.0, help="FWHM of the contraction (flux-peak strip)")
ap.add_argument("--end-time", type=float, default=4.0)
ap.add_argument("--avg-start", type=float, default=3.3)
args = ap.parse_args()

W = 0.030                                   # channel height (z), from -0.015 to +0.015
src = Path(args.src); case = Path(args.case)
times = sorted((p.name for p in src.iterdir() if re.fullmatch(r"[0-9.]+", p.name)), key=float)
t0 = times[-1]
print(f"cloning {case} from {src}/{t0}")

# flux-peak location x_peak from the fitted plasma map (argmax q over s)
q = json.loads(Path(args.qmap).read_text())
s_m = q["s_m"]; qv = q["q_W_m2"]
x_peak = s_m[max(range(len(qv)), key=lambda i: qv[i])]
f_center = 2.0 / args.speedup - 1.0
amp = 1.0 - f_center
sigma = (args.strip_mm * 1e-3) / 2.3548200            # FWHM -> Gaussian sigma
print(f"  flux peak x={x_peak*1e3:.1f} mm; target speedup {args.speedup}x -> throat height "
      f"{W*1e3*(1+f_center)/2:.2f} mm (f_center={f_center:.3f}); strip FWHM {args.strip_mm} mm")

if case.exists():
    shutil.rmtree(case)
case.mkdir()
shutil.copytree(src / "constant", case / "constant")
shutil.copytree(src / "system", case / "system")
shutil.copytree(src / t0, case / t0)

# ---- reset averaging state so the throat window is clean ----
fop = case / t0 / "uniform/functionObjects"
if fop.exists():
    shutil.rmtree(fop)
for f in set((case / t0).glob("*Mean")) | set((case / t0).glob("*Prime2Mean")):
    f.unlink()

# ---- deform constant/polyMesh/points: pinch the cold side (z<0) at the flux peak ----
pts = case / "constant/polyMesh/points"
raw = pts.read_text()
# header is everything up to and including the "(" that opens the vector list
mopen = re.search(r"\(\s*\n(?=\s*\()", raw)          # first "(" that begins the point list
head = raw[:mopen.end()]
body = raw[mopen.end():]
line_re = re.compile(r"\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)")
def deform(mobj):
    x, y, z = (float(mobj.group(i)) for i in (1, 2, 3))
    if z < 0.0:
        fx = 1.0 - amp * math.exp(-((x - x_peak) / sigma) ** 2)
        z = z * fx
    return f"({x:.10g} {y:.10g} {z:.10g})"
new_body, n = line_re.subn(deform, body)
pts.write_text(head + new_body)
print(f"  deformed {n} mesh points (cold-side z<0 pinched)")

# ---- controlDict: extend run, restart averaging after the throat re-establishes ----
cd = case / "system/controlDict"
t = cd.read_text()
t = re.sub(r"endTime\s+[0-9.]+;", f"endTime {args.end_time};", t)
t = re.sub(r"timeStart\s+[0-9.]+;", f"timeStart {args.avg_start};", t)
cd.write_text(t)

# ---- run script: NO blockMesh / NO topoSet (mesh + zones are the deformed clone) ----
rs = case / "run_throat.sh"
rs.write_text(f"""#!/usr/bin/env bash
set -e
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true
cd "$(dirname "$0")"
decomposePar -force -time {t0} > log.decomposePar 2>&1
mpirun --use-hwthread-cpus --allow-run-as-root -np 8 pimpleFoam -parallel > log.pimpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
echo DONE
""")
rs.chmod(rs.stat().st_mode | stat.S_IEXEC)

# ---- checkMesh the deformed mesh BEFORE committing to a multi-hour solve ----
chk = subprocess.run(["bash", "-c",
    f"source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1; cd {case}; "
    f"checkMesh -time {t0} 2>&1 | tail -40"], capture_output=True, text=True)
print("---- checkMesh (tail) ----")
print(chk.stdout)
