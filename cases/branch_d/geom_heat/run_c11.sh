#!/usr/bin/env bash
cd "$(dirname "$0")"
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true
decomposePar -force > log.decomposePar 2>&1
mpirun --allow-run-as-root -np 2 pimpleFoam -parallel > log.pimpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
touch C11.done
