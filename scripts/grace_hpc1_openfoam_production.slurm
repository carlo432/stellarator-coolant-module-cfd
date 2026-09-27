#!/usr/bin/env bash
#SBATCH --job-name=hpc1-of-les
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=40
#SBATCH --mem=0

set -euo pipefail

CASE_DIR="${CASE_DIR:-$PWD}"
OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-/usr/lib/openfoam/openfoam2512/etc/bashrc}"
MPI_LAUNCH="${MPI_LAUNCH:-srun}"

mkdir -p "$CASE_DIR/logs"
source "$OPENFOAM_BASHRC"
cd "$CASE_DIR"

echo "CASE_DIR=$CASE_DIR"
echo "SLURM_NTASKS=${SLURM_NTASKS:-unset}"
echo "Using OpenFOAM from $WM_PROJECT_DIR"

if [ ! -d processor0 ]; then
    decomposePar -force > log.decomposePar 2>&1
fi

$MPI_LAUNCH pimpleFoam -parallel > "log.pimpleFoam.${SLURM_JOB_ID:-manual}" 2>&1
reconstructPar -latestTime > "log.reconstructPar.${SLURM_JOB_ID:-manual}" 2>&1
postProcess -func writeCellCentres -latestTime > "log.cc.${SLURM_JOB_ID:-manual}" 2>&1 || true

echo "production chunk complete"
