#!/bin/bash
# Wall-resolved CHT velocity sweep: vary inlet U, measure converged interface wall dT.
# Tests the THERMAL side of the flushing-speed law with the trustworthy metric: dT(U) ~ U^-0.8 ?
# Reuses cht_wallres's built+split mesh; only 0/ fields (U,k,omega) change per velocity.
source /usr/lib/openfoam/openfoam2512/etc/bashrc 2>/dev/null
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RES=$W/cht_wr_sweep_results.txt
echo "# U(m/s)  Re   interfaceT  dT   Nu=(q/dT)Dh/k   heated_yplus_avg" > "$RES"
for U in 0.5799 0.8699 1.1598 1.7397 2.3196; do
  Re=$(python3 -c "print(int(round($U/1.1598*10000)))")
  K=$(python3 -c "print(round(0.00375*$U**2,6))")
  OM=$(python3 -c "import math;print(round(math.sqrt(0.00375*$U**2)/0.001023,3))")
  C=$W/cht_wr_U${U}
  rm -rf "$C"; mkdir -p "$C"
  cp -r $W/cht_wallres/system "$C/system"
  cp -r $W/cht_wallres/constant "$C/constant"      # includes the already-split polyMesh
  cp -r $W/cht_wallres/0 "$C/0"
  cd "$C"
  # set fields for this velocity
  foamDictionary -entry internalField -set "uniform ($U 0 0)" 0/bottomWater/U >/dev/null
  foamDictionary -entry 'boundaryField/inlet/value' -set "uniform ($U 0 0)" 0/bottomWater/U >/dev/null
  foamDictionary -entry internalField -set "uniform $K" 0/bottomWater/k >/dev/null
  foamDictionary -entry 'boundaryField/inlet/value' -set "uniform $K" 0/bottomWater/k >/dev/null
  foamDictionary -entry 'boundaryField/inlet/inletValue' -set "uniform $K" 0/bottomWater/k >/dev/null 2>&1
  foamDictionary -entry internalField -set "uniform $OM" 0/bottomWater/omega >/dev/null
  foamDictionary -entry 'boundaryField/inlet/value' -set "uniform $OM" 0/bottomWater/omega >/dev/null
  foamDictionary -entry 'boundaryField/inlet/inletValue' -set "uniform $OM" 0/bottomWater/omega >/dev/null 2>&1
  foamDictionary -entry stopAt -set endTime system/controlDict >/dev/null
  foamDictionary -entry endTime -set 2500 system/controlDict >/dev/null
  foamDictionary -entry writeInterval -set 2500 system/controlDict >/dev/null
  chtMultiRegionSimpleFoam > log.cht 2>&1
  IT=$(chtMultiRegionSimpleFoam -postProcess -region heater -func 'patchAverage(name=heater_to_bottomWater,fields=(T))' -latestTime 2>/dev/null | grep -oP 'of T = \K[0-9.]+')
  YP=$(chtMultiRegionSimpleFoam -postProcess -region bottomWater -func yPlus -latestTime 2>/dev/null | grep -oP 'bottomWater_to_heater y\+ : min = [0-9.]+, max = [0-9.]+, average = \K[0-9.]+')
  LINE=$(python3 -c "
dT=${IT:-800}-800
Nu=(100000/dT)*0.02667/1.0 if dT>0 else 0
print(f'$U  $Re  ${IT:-NA}  {dT:.2f}  {Nu:.1f}  ${YP:-NA}')")
  echo "$LINE" | tee -a "$RES"
done
echo "SWEEP_DONE"
echo "=== fit dT ~ U^n ==="
python3 - "$RES" <<'PY'
import sys,math
rows=[l.split() for l in open(sys.argv[1]) if l and not l.startswith('#')]
U=[float(r[0]) for r in rows]; dT=[float(r[3]) for r in rows]
# log-log linear fit dT = a*U^n
lx=[math.log(u) for u in U]; ly=[math.log(d) for d in dT]
n=len(U); sx=sum(lx); sy=sum(ly); sxx=sum(x*x for x in lx); sxy=sum(x*y for x,y in zip(lx,ly))
slope=(n*sxy-sx*sy)/(n*sxx-sx*sx)
print(f"fitted exponent n in dT ~ U^n : {slope:.3f}  (Dittus-Boelter predicts -0.8)")
print("U, dT:", list(zip(U,dT)))
PY