#!/usr/bin/env python3
"""Drift control: re-average a settled case over a fresh window on an ALREADY-DEVELOPED field.

The cold-start nowave/straight cases averaged from t=0.6 (4 flow-throughs after a uniform
initial field). Their hotspot sigma_T came out ~2x the warm-started curved case's. Either
(a) that is real -- a non-wavy duct's hotspot is set by wandering turbulent streaks rather
than pinned at a geometric weak-flushing spot -- or (b) residual thermal drift inside the
averaging window inflated TPrime2Mean (the v26 lesson: separate transient drift from
resolved fluctuation before claiming turbulence).
This clone restarts from the case's endpoint (fully developed) and averages a fresh window of
the SAME length with zero cold-start history. If sigma_hot holds, (a). If it collapses, (b).
"""
import re, shutil, stat, sys
from pathlib import Path

src, case, t_end = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3])
times = sorted((float(p.name) for p in src.iterdir() if re.fullmatch(r"[0-9.]+", p.name)), key=float)
t0s = "%g" % times[-1]
if case.exists(): shutil.rmtree(case)
case.mkdir()
for d in ("constant", "system"): shutil.copytree(src / d, case / d)
shutil.copytree(src / t0s, case / t0s)
fop = case / t0s / "uniform/functionObjects"
if fop.exists(): shutil.rmtree(fop)
for f in set((case / t0s).glob("*Mean")) | set((case / t0s).glob("*Prime2Mean")): f.unlink()
c = (case / "system/controlDict").read_text()
c = re.sub(r"endTime .*?;", f"endTime {t_end};", c, count=1)
c = re.sub(r"timeStart .*?;", f"timeStart {t0s};", c, count=1)   # average from the restart instant
(case / "system/controlDict").write_text(c)
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
rs.chmod(rs.stat().st_mode | stat.S_IEXEC)
print(f"built {case}: restart {t0s} -> {t_end}, fresh averaging window from {t0s} "
      f"({t_end - float(t0s):.2f} s = {(t_end-float(t0s))/0.15:.1f} flow-throughs), developed field")
