#!/bin/bash
# Build + run a CHT channel case at a given mesh multiplier; print interface wall T.
# Usage: bash make_cht.sh <mult> <casename>   (reuses cht_v1's working fields/materials/system)
# NOTE: source the OpenFOAM bashrc BEFORE any 'set -e' -- the bashrc returns rc=1 and
#       would otherwise kill the script at the source line (the old set-e bug).
source /usr/lib/openfoam/openfoam2512/etc/bashrc 2>/dev/null
set -u
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
M=$1; NAME=$2
C=$W/$NAME
SX=$(python3 -c "print(int(round(80*$M)))"); SY=$(python3 -c "print(int(round(8*$M)))")
NZ=$(python3 -c "print(int(round(8*$M)))"); FY=$(python3 -c "print(int(round(24*$M)))")
echo ">> building $NAME  mult=$M  solid($SX $SY $NZ)  fluid($SX $FY $NZ)"
die(){ echo "FAILED at: $1"; tail -8 "$2" 2>/dev/null; exit 1; }
rm -rf "$C"; mkdir -p "$C" || die "mkdir" /dev/null
cp -r "$W/cht_v1/system"   "$C/system"
cp -r "$W/cht_v1/constant" "$C/constant"
rm -rf "$C/constant/bottomWater/polyMesh" "$C/constant/heater/polyMesh" "$C/constant/polyMesh" 2>/dev/null
cp -r "$W/cht_v1/0" "$C/0"
rm -rf "$C/0/bottomWater/cellToRegion" "$C/0/heater/cellToRegion" 2>/dev/null
# blockMeshDict with parametrized cell counts
sed -e "s/(80  8  8)/($SX  $SY  $NZ)/" -e "s/(80 24  8)/($SX $FY  $NZ)/" \
    "$W/cht_v1/system/blockMeshDict" > "$C/system/blockMeshDict"
cd "$C" || die "cd" /dev/null
blockMesh            > log.blockMesh 2>&1 || die blockMesh log.blockMesh
grep -q "FOAM FATAL" log.blockMesh && die "blockMesh FATAL" log.blockMesh
topoSet             > log.topoSet  2>&1 || die topoSet  log.topoSet
splitMeshRegions -cellZones -overwrite > log.split 2>&1 || die split log.split
foamDictionary -entry endTime       -set 4000 system/controlDict >/dev/null
foamDictionary -entry writeInterval -set 4000 system/controlDict >/dev/null
chtMultiRegionSimpleFoam > log.cht 2>&1 || die "cht (did not converge/crashed)" log.cht
NC=$(grep -oP 'nCells:\s*\K[0-9]+' log.blockMesh | head -1)
IT=$(chtMultiRegionSimpleFoam -postProcess -region heater -func 'patchAverage(name=heater_to_bottomWater,fields=(T))' -latestTime 2>/dev/null | grep -oP 'of T = \K[0-9.]+')
BT=$(chtMultiRegionSimpleFoam -postProcess -region heater -func 'patchAverage(name=heatedBase,fields=(T))' -latestTime 2>/dev/null | grep -oP 'of T = \K[0-9.]+')
echo "CHTRESULT $NAME mult=$M cells=$NC interfaceT=${IT:-NA} baseT=${BT:-NA} (interface dT=$(python3 -c "print(round(${IT:-800}-800,3))"))"
