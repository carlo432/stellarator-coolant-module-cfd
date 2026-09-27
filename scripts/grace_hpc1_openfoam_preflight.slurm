#!/usr/bin/env bash
#SBATCH --job-name=hpc1-of-pre
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=0

set -euo pipefail

CASE_DIR="${CASE_DIR:-$PWD}"
OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-/usr/lib/openfoam/openfoam2512/etc/bashrc}"

mkdir -p "$CASE_DIR/logs"
source "$OPENFOAM_BASHRC"
cd "$CASE_DIR"

echo "CASE_DIR=$CASE_DIR"
echo "SLURM_NTASKS=${SLURM_NTASKS:-unset}"
echo "Using OpenFOAM from $WM_PROJECT_DIR"

blockMesh > log.blockMesh 2>&1
checkMesh > log.checkMesh 2>&1
setExprFields > log.setExprFields 2>&1
decomposePar -force > log.decomposePar 2>&1

echo "preflight complete"
