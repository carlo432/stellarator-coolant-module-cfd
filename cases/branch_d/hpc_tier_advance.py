#!/usr/bin/env python3
"""STEP 6 (engine) of the HPC streamlining plan: one-command tier advance.
Wires in STEP 1 (warm-start via mapFields), STEP 2 (certified recipe via gen_hpc1_channel_c4.py),
STEP 3 (WMLES vs WRLES near-wall target), STEP 6 (cost-model core sizing + SLURM emit), STEP 7 (preflight).

Generates tier N, warm-started from tier N-1's developed field, sized from the calibrated cost model, with a
Grace SLURM script and a local run script, then runs the pre-flight checker. It does NOT launch a long solver.

Example:
  # HPC-1 anchor (cold start, wall-resolved):
  python3 hpc_tier_advance.py --tier hpc1 --nx 256 --ny 192 --nz 192 --u-bar 1.16 --wrles
  # HPC-2 warm-started from hpc1, wall-modeled:
  python3 hpc_tier_advance.py --tier hpc2 --nx 320 --ny 128 --nz 256 --u-bar 1.16 --wmles --warm-start-from hpc1
"""
import argparse, math, subprocess, sys
from pathlib import Path

WORK = Path(__file__).resolve().parent
GEN = WORK / "gen_hpc1_channel_c4.py"
PREFLIGHT = WORK / "hpc_preflight_check.py"
COST_C = 1.06e-5          # core-s = C * Ncells^(4/3) * flowthroughs  (calibrated locally)

def divisors(n):
    return [d for d in range(1, n+1) if n % d == 0]

def pick_decomp(nx, ny, nz, target_cpc):
    """cost-model core sizing: aim ~target_cpc cells/core; split x,z (homogeneous), minimal y."""
    ncells = nx*ny*nz
    target_cores = max(1, round(ncells/target_cpc))
    best = None
    for dy in (1, 2, 4):
        if ny % dy: continue
        for dx in divisors(nx):
            for dz in divisors(nz):
                cores = dx*dy*dz
                if cores < 1: continue
                cpc = ncells/cores
                if not (0.5*target_cpc <= cpc <= 2.0*target_cpc): continue
                # prefer cores near target, dx~dz (balanced homogeneous split), small dy
                score = (abs(cores-target_cores), abs(dx-dz), dy)
                if best is None or score < best[0]:
                    best = (score, (dx, dy, dz), cores, cpc)
    if best is None:  # fallback: dy=1, dx=dz nearest
        c = max(1, round(math.sqrt(target_cores)))
        return (c, 1, max(1, target_cores//c)), c*max(1, target_cores//c), ncells/(c*max(1,target_cores//c))
    return best[1], best[2], best[3]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", required=True)
    ap.add_argument("--nx", type=int, required=True); ap.add_argument("--ny", type=int, required=True); ap.add_argument("--nz", type=int, required=True)
    ap.add_argument("--u-bar", type=float, default=1.1598); ap.add_argument("--nu", type=float, default=3.09278e-6); ap.add_argument("--pr", type=float, default=14.4)
    ap.add_argument("--wmles", action="store_true", help="wall-modeled (coarser near-wall grading)")
    ap.add_argument("--wrles", action="store_true", help="wall-resolved (fine near-wall grading)")
    ap.add_argument("--warm-start-from", default=None, help="STEP 1: source tier case to mapFields developed field from")
    ap.add_argument("--flowthroughs", type=float, default=20.0, help="for the cost-model core-hour estimate")
    ap.add_argument("--cells-per-core", type=float, default=35000.0)
    ap.add_argument("--no-thermal", action="store_true")
    args = ap.parse_args()

    case = f"tier_{args.tier}"
    ncells = args.nx*args.ny*args.nz
    grad_y = 14.0 if args.wmles else 36.0     # STEP 3: WMLES coarser near-wall; WRLES fine (y+<1, thermal-aware)
    (dx, dy, dz), cores, cpc = pick_decomp(args.nx, args.ny, args.nz, args.cells_per_core)
    core_s = COST_C * ncells**(4.0/3.0) * args.flowthroughs
    core_hr = core_s/3600.0

    # STEP 2: generate with the certified recipe (gen carries WALE + alphaD/alphaDt thermal + gotcha-safe dicts)
    cmd = [sys.executable, str(GEN), "--case", case, "--nx", str(args.nx), "--ny", str(args.ny), "--nz", str(args.nz),
           "--cores", str(cores), "--decomp", str(dx), str(dy), str(dz), "--u-bar", str(args.u_bar),
           "--nu", str(args.nu), "--pr", str(args.pr), "--grad-y", str(grad_y)]
    if args.no_thermal: cmd.append("--no-thermal")
    print("generating:", " ".join(cmd))
    subprocess.run(cmd, cwd=WORK, check=True)

    # STEP 1: warm-start script (mapFields the developed field of the source tier)
    ws = ""
    if args.warm_start_from:
        src = f"tier_{args.warm_start_from}" if not args.warm_start_from.startswith("tier_") else args.warm_start_from
        ws = (f"# STEP 1 warm-start: interpolate developed field from {src} (skip transition).\n"
              f"# (interpolating map -- tiers have DIFFERENT meshes; -consistent is only for identical meshes)\n"
              f"mapFields ../{src} -case . -sourceTime latestTime > log.mapFields 2>&1\n")

    # local run script (gotcha #7) and Grace SLURM (srun)
    (WORK/case/"run_local.sh").write_text(
        "#!/usr/bin/env bash\nexport OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1\n"
        "source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true\ncd \"$(dirname \"$0\")\"\n"
        "blockMesh > log.blockMesh 2>&1\nsetExprFields > log.setExprFields 2>&1\ncheckMesh > log.checkMesh 2>&1\n"
        + ws +
        f"decomposePar -force > log.decomposePar 2>&1\n"
        f"mpirun --allow-run-as-root -np {cores} pimpleFoam -parallel > log.pimpleFoam 2>&1\n"
        "reconstructPar -latestTime > log.reconstructPar 2>&1\n")
    (WORK/case/"submit_grace.slurm").write_text(
        f"#!/bin/bash\n#SBATCH -J {case}\n#SBATCH -N {max(1, math.ceil(cores/48))}\n#SBATCH -n {cores}\n"
        f"#SBATCH -t 24:00:00\n#SBATCH -p cpu\nmodule load OpenFOAM/v2512 || true\nsource $FOAM_BASHRC\n"
        "blockMesh > log.blockMesh 2>&1\nsetExprFields > log.setExprFields 2>&1\n"
        + (ws.replace("mapFields", "mapFields") if ws else "") +
        f"decomposePar -force > log.decomposePar 2>&1\nsrun -n {cores} pimpleFoam -parallel > log.pimpleFoam 2>&1\n"
        "reconstructPar -latestTime > log.reconstructPar 2>&1\n")

    print("\n=== TIER PLAN ===")
    print(f"  tier         : {args.tier}  ({'WMLES' if args.wmles else 'WRLES'}, grad_y={grad_y})")
    print(f"  cells        : {ncells:,} ({args.nx}x{args.ny}x{args.nz})")
    print(f"  decomp       : {cores} cores  n=({dx} {dy} {dz})  -> {cpc:,.0f} cells/core")
    print(f"  warm-start   : {'from '+args.warm_start_from+' (mapFields)' if args.warm_start_from else 'COLD (no source)'}")
    print(f"  cost model   : ~{core_hr:,.0f} core-hr  ({args.flowthroughs:.0f} flow-throughs)")
    print(f"  scripts      : {case}/run_local.sh , {case}/submit_grace.slurm")

if __name__ == "__main__":
    main()
