#!/usr/bin/env python3
"""Warm-start the straight-pipeline PLASMA leg: clone geom_straight_pipeline at its latest
time, swap the heated-wall exprMixed to the plasma-shaped map with s = pos().x().

The curved plasma map fits the smooth standoff d(s) as a polynomial in t = 2s/L - 1, with
s = atan2(x, 0.18-z)*0.18 (arc-length unroll). For the straight duct L is identical and
s = pos().x() directly, so the SAME fitted polynomial applies with that one substitution --
identical delivered load shape, 1:1, verifying the amplification isn't a curved-unroll
artifact either.
"""
from __future__ import annotations
import argparse, json, re, shutil
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--src", default="geom_straight_pipeline")
ap.add_argument("--case", default="geom_straight_pipeline_plasma")
ap.add_argument("--qmap", default="figs/plasma_qmap.json")
ap.add_argument("--end-time", type=float, default=2.7)
ap.add_argument("--avg-start", type=float, default=2.0)
ap.add_argument("--coord", choices=("x", "arc"), default="x",
                help="streamwise coordinate for the map: 'x' (straight duct) or 'arc' "
                     "(curved duct -- keep the atan2 unroll verbatim)")
args = ap.parse_args()

src = Path(args.src)
case = Path(args.case)
# latest reconstructed time dir
times = sorted((float(p.name) for p in src.iterdir() if re.fullmatch(r"[0-9.]+", p.name)),
               key=float)
t0 = times[-1]
t0s = ("%g" % t0)
print(f"warm-starting {case} from {src}/{t0s}")

if case.exists():
    shutil.rmtree(case)
case.mkdir()
shutil.copytree(src / "constant", case / "constant")
shutil.copytree(src / "system", case / "system")
shutil.copytree(src / t0s, case / t0s)
# reset averaging state, drop old Mean/Prime2Mean so the plasma window is clean
fop = case / t0s / "uniform/functionObjects"
if fop.exists():
    shutil.rmtree(fop)
for f in set((case / t0s).glob("*Mean")) | set((case / t0s).glob("*Prime2Mean")):
    f.unlink()

# straight plasma gradientExpr: curved arc-unroll -> pos().x()
q = json.loads(Path(args.qmap).read_text())
gexpr = q["gradientExpr"]
if args.coord == "x":
    gexpr = gexpr.replace("atan2(pos().x(), 0.18 - pos().z()) * 0.18", "pos().x()")
    assert "atan2" not in gexpr, "arc-unroll substitution incomplete"
else:
    assert "atan2" in gexpr, "curved map should retain the arc unroll"

tfile = case / t0s / "T"
txt = tfile.read_text()
old = re.search(r"    heated_first_wall\n    \{.*?\n    \}\n", txt, re.S).group(0)
new = ("    heated_first_wall\n    {\n"
       "        type            exprMixed;\n"
       "        refValue        uniform 900;\n"
       "        refGradient     uniform 500646;\n"
       "        valueFraction   uniform 0;\n"
       "        source          uniform 0;\n"
       "        value           uniform 900;\n"
       "        valueExpr       #{ 900 #};\n"
       f"        gradientExpr    #{{ {gexpr} #}};\n"
       "        fractionExpr    #{ 0 #};\n"
       "        variables       ();\n"
       "    }\n")
tfile.write_text(txt.replace(old, new))

# controlDict: extend endTime + reset averaging window start
cd = case / "system/controlDict"
c = cd.read_text()
c = re.sub(r"endTime .*?;", f"endTime {args.end_time};", c, count=1)
c = re.sub(r"timeStart .*?;", f"timeStart {args.avg_start};", c, count=1)
cd.write_text(c)

rs = case / "run_pipeline.sh"
rs.write_text(f"""#!/usr/bin/env bash
set -e
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true
cd "$(dirname "$0")"
decomposePar -force -time {t0s} > log.decomposePar 2>&1
mpirun --use-hwthread-cpus --allow-run-as-root -np 8 pimpleFoam -parallel > log.pimpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
echo DONE
""")
import stat as _s
rs.chmod(rs.stat().st_mode | _s.S_IEXEC)
coord_desc = "s=pos().x()" if args.coord == "x" else "s=arc unroll (atan2)"
print(f"  plasma map applied ({coord_desc}); endTime {args.end_time}, avg from {args.avg_start}")
qpk = q.get("q_peak_W_m2")
extra = f", q_peak {qpk/1e6:.3f} MW/m2" if qpk else ""
print(f"  peak/mean {q['peak_over_mean']:.3f}{extra}  [map: {args.qmap}]")
