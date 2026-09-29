#!/usr/bin/env python3
"""THE FUSE -- one-command GEOMETRY tier (the non-periodic analog of hpc_tier_advance.py).

Channel tiers (HPC-1/2, periodic) use hpc_tier_advance.py. Geometry tiers (HPC-3..6: curved slice,
multi-segment, stellarator-CAD) are NON-periodic, so they need a different spine. This fuses the four
geometry-relevant streamlining steps into one driver:

  base mesh  : gen_curved_stellarator_slice_c2.py            (the curved/stellarator geometry + patches)
  LES+thermal: overlay RAS-smoke -> WALE LES + transported T (alphaD/alphaDt SGS flux, the C4 fix)
  STEP 4     : hpc_inflow_dfsem.py   precursor -> turbulentDFSEM inlet (no development length grown in-domain)
  STEP 5     : hpc_freeze_neutronics.py  OpenMC q''' -> frozen scalarSemiImplicitSource on T (run OpenMC once)
  STEP 6     : cost-model core sizing + run_local.sh / submit_grace.slurm
  STEP 7     : hpc_preflight_check.py  (catch the 7 gotchas before the queue)
  STEP 1     : optional warm-start via mapFields from a prior geometry tier

This is exactly the C11 synthesis (curved slice + C1 thermal recipe), produced as the fuse's first output.

NOTE on the two honest scaling knobs (the fuse wires the mechanics; these still need a human eye):
  * inlet remap: the precursor's half-height delta should match the slice gap; DFSEM nearestCell will clamp
    otherwise. Pass --precursor-delta or accept the precursor's own delta.
  * q''' -> K/s: a transported T [K] needs a source in K/s = q'''/(rho*cp). --rho-cp converts; the C3 CSV
    'normalized_qvol' is a SHAPE (per-source-particle), so --q-scale sets the absolute peak W/m3.

Usage:
  python3 hpc_geometry_tier.py --tier c11 --precursor les_channel_c8_hiRe \
      --neutronics figs/c3_openmc_flibe_slab_heating_profile.csv --q-scale 1.5e7 \
      --segments 24 --radial-cells 24 --height-cells 24 --stream-cells 64 \
      --y-amp 0.0

Geometry note:
  --y-amp is intentionally explicit.  The earlier 10 mm sinusoidal span waviness
  costs pumping and raises the uniform-load mean hotspot, but drift-controlled
  flicker analysis showed it can halve plasma-hotspot sigma_T.  Mean-film-
  temperature-limited studies may prefer y_amp=0; transient/fatigue-limited
  studies may keep waviness.  Do not hide that trade behind a default.
"""
import argparse, math, os, re, shutil, subprocess, sys
from pathlib import Path

WORK = Path(__file__).resolve().parent
GEN = WORK / "gen_curved_stellarator_slice_c2.py"
DFSEM = WORK / "hpc_inflow_dfsem.py"
FREEZE = WORK / "hpc_freeze_neutronics.py"
PREFLIGHT = WORK / "hpc_preflight_check.py"
# geometry-tier cost coefficient, CALIBRATED on C11 (curved + thermal scalar + DFSEM): ~18.5x the periodic-channel
# C=1.06e-5. Using the channel value here under-provisions cores ~18x -> jobs blow the Grace allocation.
COST_C = 1.96e-4

def kader_beta(pr):
    # Kader (1981): beta(Pr) = (3.85*Pr^(1/3) - 1.3)^2 + 2.12*ln(Pr)
    return (3.85 * pr ** (1.0 / 3.0) - 1.3) ** 2 + 2.12 * math.log(pr)

def thermal_wall_boundary(args):
    if args.thermal_wallmodel == "none":
        return '    "(heated_first_wall|blanket_back_wall|side_walls)" { type zeroGradient; }\n'
    if args.thermal_wallmodel == "kader-static":
        theta_plus = 2.12 * math.log(args.thermal_wallmodel_yplus) + kader_beta(args.pr)
        wall_t = args.t_in + args.thermal_wallmodel_theta_scale * theta_plus
        return (
            '    "(heated_first_wall|blanket_back_wall|side_walls)"\n'
            "    {\n"
            "        type fixedValue;\n"
            f"        // OPT-1 static Kader surrogate: Pr={args.pr:.8g}, yPlusRef={args.thermal_wallmodel_yplus:.8g}, thetaPlus={theta_plus:.8g}\n"
            f"        // Smoke/prototype only; not the final mixed thermal wall function.\n"
            f"        value uniform {wall_t:.12g};\n"
            "    }\n"
        )
    raise ValueError(f"unknown thermal wall model: {args.thermal_wallmodel}")

def divisors(n): return [d for d in range(1, n + 1) if n % d == 0]

def pick_decomp(ncells, nx, ny, nz, target_cpc):
    """split stream(nx) x height(nz); minimal radial(ny) -- homogeneous-ish directions get the cuts."""
    target_cores = max(1, round(ncells / target_cpc))
    best = None
    for dy in (1, 2):
        if ny % dy: continue
        for dx in divisors(nx):
            for dz in divisors(nz):
                cores = dx * dy * dz
                cpc = ncells / cores
                if not (0.5 * target_cpc <= cpc <= 2.0 * target_cpc): continue
                score = (abs(cores - target_cores), abs(dx - dz), dy)
                if best is None or score < best[0]:
                    best = (score, (dx, dy, dz), cores, cpc)
    if best is None:
        c = max(1, round(math.sqrt(target_cores)))
        d = max(1, target_cores // c)
        return (c, 1, d), c * d, ncells / (c * d)
    return best[1], best[2], best[3]

def header(cls, obj):
    return ("FoamFile\n{\n    version 2.0; format ascii; class %s; object %s;\n}\n" % (cls, obj))

def overlay_les_thermal(case, args, freeze_block):
    """rewrite the RAS-smoke dicts into a WALE LES with a transported, neutronically-heated scalar T."""
    c = Path(case)
    # constant/turbulenceProperties -> LES WALE
    (c / "constant/turbulenceProperties").write_text(
        header("dictionary", "turbulenceProperties") +
        "simulationType LES;\nLES\n{\n    LESModel WALE;\n    turbulence on;\n    printCoeffs on;\n"
        "    delta cubeRootVol;\n    cubeRootVolCoeffs { deltaCoeff 1; }\n}\n")
    # constant/transportProperties keep nu; DT unused (scalarTransport uses alphaD/alphaDt path)
    # 0/T : real inlet/outlet temperature scalar (K), walls adiabatic (volumetric nuclear heat)
    wall_bc = thermal_wall_boundary(args)
    (c / "0/T").write_text(
        header("volScalarField", "T") +
        "dimensions [0 0 0 1 0 0 0];\n"
        f"internalField uniform {args.t_in:.6g};\n"
        "boundaryField\n{\n"
        f"    inlet {{ type fixedValue; value uniform {args.t_in:.6g}; }}\n"
        f"    outlet {{ type inletOutlet; inletValue uniform {args.t_in:.6g}; value uniform {args.t_in:.6g}; }}\n"
        f"{wall_bc}"
        "}\n")
    # system/controlDict -> pimpleFoam, physical time, fieldAverage(writeTime), scalarTransport(T)+frozen source
    nu = args.nu; Pr = args.pr; Prt = 0.9
    alphaD = 1.0 / Pr; alphaDt = 1.0 / Prt
    (c / "system/controlDict").write_text(
        header("dictionary", "controlDict") +
        f"""application pimpleFoam;
startFrom latestTime;
startTime 0;
stopAt endTime;
endTime {args.end_time:.6g};
deltaT {args.dt:.6g};
writeControl adjustableRunTime;
writeInterval {args.write_interval:.6g};
purgeWrite 3;
writeFormat ascii;
writePrecision 8;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable true;
adjustTimeStep yes;
maxCo 0.6;
maxDeltaT {args.dt*5:.6g};
functions
{{
    fieldAverage1
    {{
        type fieldAverage; libs (fieldFunctionObjects);
        writeControl writeTime;
        timeStart {args.avg_start:.6g};
        fields ( U {{ mean on; prime2Mean on; base time; }} T {{ mean on; prime2Mean on; base time; }} );
    }}
    transportT
    {{
        type scalarTransport; libs (solverFunctionObjects);
        field T; schemesField U; resetOnStartUp false;
        // alphaD/alphaDt force the turbulence-model path: D_eff = alphaD*nu + alphaDt*nut
        alphaD {alphaD:.8g}; alphaDt {alphaDt:.8g};
        fvOptions
        {{
{freeze_block}
        }}
    }}
    yPlus1 {{ type yPlus; libs (fieldFunctionObjects); writeControl writeTime; }}
}}
""")
    # system/fvSchemes -> LES (backward ddt, LUST momentum, bounded scalar)
    (c / "system/fvSchemes").write_text(
        header("dictionary", "fvSchemes") +
        """ddtSchemes { default backward; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) Gauss LUST grad(U);
    div(phi,T) Gauss limitedLinear 1;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    # system/fvSolution -> PIMPLE + UFinal + TFinal (no pRefCell: outlet fixes pressure)
    (c / "system/fvSolution").write_text(
        header("dictionary", "fvSolution") +
        """solvers
{
    p { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0.01; }
    pFinal { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0; }
    "(U|T)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0.1; }
    "(U|T)Final" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0; }
}
PIMPLE
{
    nOuterCorrectors 1;
    nCorrectors 2;
    nNonOrthogonalCorrectors """ + str(args.non_ortho_correctors) + """;
}
""")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="c11")
    ap.add_argument("--precursor", default=None, help="STEP 4: developed periodic-channel case for DFSEM inflow")
    ap.add_argument("--precursor-delta", default="auto")
    ap.add_argument("--no-inflow", action="store_true", help="skip DFSEM; keep laminar fixedValue inlet")
    ap.add_argument("--neutronics", default=None, help="STEP 5: OpenMC heating CSV (C3 format)")
    ap.add_argument("--q-scale", type=float, default=None, help="absolute peak q''' [W/m3] to scale the C3 shape to")
    ap.add_argument("--rho-cp", type=float, default=4.65e6, help="FLiBe rho*cp [J/m3/K] for q'''->K/s")
    ap.add_argument("--warm-start-from", default=None, help="STEP 1: prior geometry tier to mapFields from")
    # geometry (passed through to the curved generator)
    ap.add_argument("--segments", type=int, default=24); ap.add_argument("--stream-cells", type=int, default=64)
    ap.add_argument("--radial-cells", type=int, default=24); ap.add_argument("--height-cells", type=int, default=24)
    ap.add_argument("--radius", type=float, default=0.18); ap.add_argument("--arc-deg", type=float, default=95.0)
    ap.add_argument("--y-amp", type=float, required=True,
                    help="span waviness amplitude [m]; explicit design choice, e.g. 0.0 for clean mean-T baseline or 0.01 for the prior wavy/flicker-damped geometry")
    # physics / run
    ap.add_argument("--u-in", type=float, default=0.86); ap.add_argument("--nu", type=float, default=3.09278e-6)
    ap.add_argument("--pr", type=float, default=14.4); ap.add_argument("--t-in", type=float, default=900.0)
    ap.add_argument("--end-time", type=float, default=0.6); ap.add_argument("--dt", type=float, default=2e-4)
    ap.add_argument("--write-interval", type=float, default=0.05)
    ap.add_argument("--avg-start", type=float, default=None,
                    help="fieldAverage timeStart; default = 3 flow-throughs (ramp ~2FT + 1FT margin) from "
                         "arc length/u_in. NOTE: u_in is the LABEL; the banked C8 inlet delivers ~1.99 m/s "
                         "(OPERATING_POINT.md), so the label-based default is conservative (starts later).")
    ap.add_argument("--flowthroughs", type=float, default=15.0); ap.add_argument("--cells-per-core", type=float, default=35000.0)
    ap.add_argument("--non-ortho-correctors", type=int, default=1,
                    help="OPT-2 profile knob: PIMPLE nNonOrthogonalCorrectors, test 0 vs 1 on low-skew geometry")
    ap.add_argument("--thermal-wallmodel", choices=("none", "kader-static"), default="none",
                    help="OPT-1 prototype knob; 'kader-static' emits a fixedValue Kader surrogate for smoke tests only")
    ap.add_argument("--thermal-wallmodel-yplus", type=float, default=30.0,
                    help="reference y+ for the static Kader surrogate")
    ap.add_argument("--thermal-wallmodel-theta-scale", type=float, default=0.1,
                    help="K per theta+ for the static Kader surrogate")
    ap.add_argument("--smoke-steps", type=int, default=8, help="GATE A: serial runnability smoke step count")
    ap.add_argument("--no-smoke", action="store_true", help="skip the runnability smoke gate")
    args = ap.parse_args()
    args.dt = args.dt  # keep
    if args.avg_start is None:
        ft = math.radians(args.arc_deg) * args.radius / args.u_in
        args.avg_start = 3.0 * ft
        print(f"avg-start default: 3 flow-throughs = {args.avg_start:.4g} s (FT={ft:.4g} s at label u_in)")
    case = WORK / f"geom_{args.tier}"

    # --- base curved geometry (the mesh + patches), RAS-smoke we then overlay ---
    gcmd = [sys.executable, str(GEN), "--case", str(case), "--segments", str(args.segments),
            "--stream-cells", str(args.stream_cells), "--radial-cells", str(args.radial_cells),
            "--height-cells", str(args.height_cells), "--radius", str(args.radius),
            "--y-amp", str(args.y_amp),
            "--arc-deg", str(args.arc_deg), "--u-in", str(args.u_in), "--nu", str(args.nu)]
    print("base geometry:", " ".join(gcmd))
    subprocess.run(gcmd, check=True)

    # along-stream cells = segments*stream-cells (the generator sweeps segments)
    nx = args.segments * args.stream_cells; ny = args.radial_cells; nz = args.height_cells
    ncells = nx * ny * nz

    # --- STEP 5: freeze neutronics -> a scalarSemiImplicitSource block (indented for the FO fvOptions) ---
    freeze_block = ("            // STEP 5 skipped (no --neutronics): inactive zero source kept syntactically complete\n"
                    "            none\n"
                    "            {\n"
                    "                type scalarSemiImplicitSource;\n"
                    "                active false;\n"
                    "                selectionMode all;\n"
                    "                volumeMode specific;\n"
                    "                sources { T (0 0); }\n"
                    "            }")
    if args.neutronics:
        # convert q'''[W/m3] -> source on T [K/s] via rho*cp; scale C3 shape to --q-scale peak if given
        fcmd = [sys.executable, str(FREEZE), str(WORK / args.neutronics) if not Path(args.neutronics).is_absolute()
                else args.neutronics, "--field", "T", "--mode", "uniform"]
        fr = subprocess.run(fcmd, cwd=WORK, capture_output=True, text=True)
        print(fr.stdout.strip())
        # peak/mean ~1 for the FLiBe slab -> uniform K/s source = q_scale/(rho*cp) (or CSV mean if no scale)
        import csv as _csv
        rows = list(_csv.DictReader(open(WORK / args.neutronics if not Path(args.neutronics).is_absolute() else args.neutronics)))
        qmean = sum(float(r["normalized_qvol_W_m3"]) for r in rows) / len(rows)
        qpeak = max(float(r["normalized_qvol_W_m3"]) for r in rows)
        qabs = args.q_scale if args.q_scale else qmean        # absolute W/m3
        src_K_s = (qabs / qpeak) * qmean / args.rho_cp if args.q_scale else qmean / args.rho_cp
        if args.q_scale:  # scale shape so peak == q_scale, use mean for the (near-flat) uniform source
            src_K_s = (args.q_scale * (qmean / qpeak)) / args.rho_cp
        freeze_block = (
            "            frozenNeutronicsHeating\n            {\n"
            "                type scalarSemiImplicitSource; active true; selectionMode all;\n"
            "                volumeMode specific;   // K/s = q'''/(rho*cp), frozen OpenMC shape (C3 near-flat)\n"
            f"                sources {{ T ({src_K_s:.6g} 0); }}\n            }}")
        print(f"  STEP5: q_abs={qabs:.4g} W/m3, rho*cp={args.rho_cp:.3g} -> T source {src_K_s:.4g} K/s")

    # --- overlay LES + thermal onto the base case ---
    overlay_les_thermal(case, args, freeze_block)

    # --- STEP 4: DFSEM inflow. A geometry inlet is a 2D PLANE, so blockMesh first, then dfsem_remap_2d.py
    # builds boundaryData at the real face centres (1D-line boundaryData fails -- 'points on a single line'),
    # and we emit the VERIFIED turbulentDFSEMInlet BC (U/R/L as PatchFunction1 mappedFile, from feature_smoke). ---
    inflow_note = "inlet: laminar fixedValue (no --precursor / --no-inflow)"
    if args.precursor and not args.no_inflow:
        pre = WORK / args.precursor if not Path(args.precursor).is_absolute() else Path(args.precursor)
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        subprocess.run(["bash","-lc","source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true; "
                        "blockMesh > log.blockMesh 2>&1"], cwd=case, env=env)
        rr = subprocess.run([sys.executable, str(WORK/"dfsem_remap_2d.py"), str(case), str(pre), "--patch", "inlet"],
                            capture_output=True, text=True)
        print(rr.stdout.strip() or rr.stderr.strip())
        dm = re.search(r"delta_chan=([0-9.eE+-]+)", rr.stdout)
        delta = float(dm.group(1)) if dm else 0.008
        if (case/"constant"/"boundaryData"/"inlet"/"points").exists():
            uf = case / "0" / "U"; txt = uf.read_text()
            bc = (f"inlet\n    {{\n        type            turbulentDFSEMInlet;\n        delta           {delta:.6g};\n"
                  "        nCellPerEddy    5;\n        U               { type mappedFile; }\n"
                  "        R               { type mappedFile; }\n        L               { type mappedFile; }\n"
                  f"        value           uniform ({args.u_in:.6g} 0 0);\n    }}")
            txt = re.sub(r"inlet\s*\{[^}]*\}", bc, txt, count=1)
            uf.write_text(txt)
            inflow_note = f"inlet: turbulentDFSEMInlet (2D-remapped plane) from {pre.name}, delta={delta:.4g}"

    # --- STEP 1: warm-start script line ---
    ws = ""
    if args.warm_start_from:
        src = f"geom_{args.warm_start_from}" if not args.warm_start_from.startswith("geom_") else args.warm_start_from
        ws = (f"# STEP 1 warm-start (interpolating; different mesh): \n"
              f"mapFields ../{src} -case . -sourceTime latestTime > log.mapFields 2>&1\n")

    # --- STEP 6: cost-model decomp + run scripts ---
    (dx, dy, dz), cores, cpc = pick_decomp(ncells, nx, ny, nz, args.cells_per_core)
    core_hr = COST_C * ncells ** (4.0 / 3.0) * args.flowthroughs / 3600.0
    # OPT-3: curved geometry meshes do NOT fill their bounding box uniformly, so geometric `simple`
    # decomposition load-imbalances badly (measured 3.92x: slowest core 144k cells vs 36.8k avg on the
    # 884k slice) and moves 5.5x more MPI traffic. `scotch` graph-partitions to perfect cell balance and
    # minimal edge cut -- ~3.9x faster wall-clock at no accuracy cost, and it is the right default for the
    # curved tier. `simple` is fine ONLY for the box channel precursors. Both are local + Grace available.
    (case / "system" / "decomposeParDict").write_text(
        header("dictionary", "decomposeParDict") +
        f"numberOfSubdomains {cores};\nmethod scotch;\n"
        f"// simple fallback (box meshes only): method simple; coeffs {{ n ({dx} {dy} {dz}); }}\n")
    (case / "run_local.sh").write_text(
        "#!/usr/bin/env bash\nexport OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1\n"
        "source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true\ncd \"$(dirname \"$0\")\"\n"
        "blockMesh > log.blockMesh 2>&1\ncheckMesh > log.checkMesh 2>&1\n" + ws +
        f"decomposePar -force > log.decomposePar 2>&1\n"
        f"mpirun --allow-run-as-root -np {cores} pimpleFoam -parallel > log.pimpleFoam 2>&1\n"
        "reconstructPar -latestTime > log.reconstructPar 2>&1\n")
    (case / "submit_grace.slurm").write_text(
        f"#!/bin/bash\n#SBATCH -J geom_{args.tier}\n#SBATCH -N {max(1, math.ceil(cores/48))}\n#SBATCH -n {cores}\n"
        f"#SBATCH -t 24:00:00\n#SBATCH -p cpu\nmodule load OpenFOAM/v2512 || true\nsource $FOAM_BASHRC\n"
        "blockMesh > log.blockMesh 2>&1\n" + ws +
        f"decomposePar -force > log.decomposePar 2>&1\nsrun -n {cores} pimpleFoam -parallel > log.pimpleFoam 2>&1\n"
        "reconstructPar -latestTime > log.reconstructPar 2>&1\n")

    print("\n=== GEOMETRY TIER (the fuse) ===")
    print(f"  tier        : {args.tier}  (curved slice, WALE LES + transported T)")
    print(f"  cells       : {ncells:,} ({nx} stream x {ny} radial x {nz} height)")
    print(f"  geometry    : radius={args.radius:g} m, arc={args.arc_deg:g} deg, y_amp={args.y_amp:g} m")
    print(f"  decomp      : {cores} cores n=({dx} {dy} {dz}) -> {cpc:,.0f} cells/core")
    print(f"  STEP4 {inflow_note}")
    print(f"  STEP5 {'frozen neutronics on T' if args.neutronics else 'no neutronics'}")
    print(f"  STEP1 warm-start: {'from '+args.warm_start_from if args.warm_start_from else 'COLD'}")
    print(f"  PIMPLE     : nNonOrthogonalCorrectors={args.non_ortho_correctors}")
    print(f"  thermal WF : {args.thermal_wallmodel}")
    print(f"  cost model  : ~{core_hr:,.0f} core-hr ({args.flowthroughs:.0f} flow-throughs)")
    print(f"  scripts     : {case.name}/run_local.sh , submit_grace.slurm")

    # --- STEP 7: preflight (structural) ---
    print()
    subprocess.run([sys.executable, str(PREFLIGHT), str(case)], check=True)

    # --- GATE B: inlet contract (Re match / 1D-vs-2D / delta<->gap) if a precursor was wired ---
    if args.precursor and not args.no_inflow:
        print()
        subprocess.run([sys.executable, str(WORK/"inlet_contract.py"), str(case),
                        str(WORK/args.precursor if not Path(args.precursor).is_absolute() else args.precursor),
                        "--u-tier", str(args.u_in), "--nu", str(args.nu)], check=True)

    # --- GATE A: short serial RUNNABILITY smoke (the real catch -- preflight proves shape, not that it runs) ---
    if not args.no_smoke:
        print(f"\n=== runnability smoke ({args.smoke_steps} steps, serial) ===")
        cd = case/"system"/"controlDict"; orig = cd.read_text()
        smk = re.sub(r"endTime\s+[0-9.eE+-]+;", f"endTime {args.smoke_steps*args.dt:.8g};", orig)
        smk = re.sub(r"startFrom\s+\w+;", "startFrom startTime;", smk)
        smk = re.sub(r"adjustTimeStep\s+\w+;", "adjustTimeStep no;", smk)
        cd.write_text(smk)
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        bcmd = (f"source /usr/lib/openfoam/openfoam2512/etc/bashrc >/dev/null 2>&1 || true; "
                f"blockMesh > log.smoke 2>&1; pimpleFoam >> log.smoke 2>&1")
        try:
            subprocess.run(["bash","-lc",bcmd], cwd=case, env=env, timeout=300)
        except subprocess.TimeoutExpired:
            pass
        cd.write_text(orig)                                    # restore the real controlDict
        log = (case/"log.smoke").read_text() if (case/"log.smoke").exists() else ""
        fatal = "FOAM FATAL" in log
        adv = len(re.findall(r"^Time = ", log, re.M)) >= 2
        tsolve = "Solving for T" in log
        if fatal:
            em = re.search(r"-->\s*FOAM FATAL.*?\n(.*?)(?:\n\n|\nFrom )", log, re.S)
            print(f"  [FAIL] solver crashed: {(em.group(1).strip()[:140] if em else '?')}")
            print(f"         -> see {case.name}/log.smoke  (DO NOT submit to Grace until this is green)")
        elif adv and tsolve:
            print(f"  [ OK ] solver advanced {len(re.findall(r'^Time = ',log,re.M))} steps, T transported -> case is RUNNABLE")
        else:
            print(f"  [WARN] no clean time advance / T not solved -> inspect {case.name}/log.smoke")

if __name__ == "__main__":
    main()
