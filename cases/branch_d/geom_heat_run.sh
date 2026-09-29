#!/bin/bash
# HEAT FLOW for the stellarator -- the thing Leffler listed as future work, on our 3D curved slice.
# Warm-start from C11's developed flow (t=0.6), apply a realistic 0.5 MW/m2 plasma heat flux on the
# heated first wall, develop the thermal field. scotch-decomposed, parallel.
set -e
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
cd "$(dirname "${BASH_SOURCE[0]}")"
NP=8
SN=500650           # snGrad(T) for q''=0.5 MW/m2 (rho_cp*D_eff ~ 1)
rm -rf geom_heat; cp -r geom_c11 geom_heat
cd geom_heat
rm -rf processor*

# 1) heated first-wall BC: zeroGradient -> fixedGradient (plasma load) in the restart time dir 0.6/T
python3 - "$SN" <<'PY'
import re, sys
sn = sys.argv[1]
f = "0.6/T"; s = open(f).read()
# target the explicit heated_first_wall block (solver writes each wall separately)
pat = re.compile(r'(heated_first_wall\s*\{\s*)type\s+zeroGradient;(\s*\})', re.S)
repl = r'\1type            fixedGradient;\n        gradient        uniform ' + sn + r';\2'
s2, n = pat.subn(repl, s)
assert n == 1, f"expected 1 heated_first_wall match, got {n}"
open(f, "w").write(s2)
print("  heated_first_wall -> fixedGradient", sn)
PY

# 2) continue from t=0.6, develop thermal field to t=0.8
sed -i 's/^startFrom .*/startFrom latestTime;/' system/controlDict
sed -i 's/^endTime .*/endTime 0.8;/' system/controlDict
sed -i 's/^writeInterval .*/writeInterval 0.05;/' system/controlDict

# 3) scotch decomposition (OPT-3) on NP cores
cat > system/decomposeParDict <<EOF
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains ${NP};
method scotch;
EOF

decomposePar -latestTime -force > log.decomposePar 2>&1
mpirun --use-hwthread-cpus --allow-run-as-root -np ${NP} pimpleFoam -parallel > log.pimpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
touch GEOM_HEAT.done
echo "geom_heat thermal run complete -> $(ls -d [0-9]* | sort -g | tail -1)"
