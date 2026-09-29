#!/usr/bin/env python3
"""STEP 7 of the HPC streamlining plan: pre-flight a case against the 7 pipeline gotchas BEFORE it touches the
Grace queue. Read-only; scans system/ dicts + checkMesh log. Exit 0 if no FAILs.

Usage: hpc_preflight_check.py <case_dir>
"""
import json, math, re, sys, os, glob
from pathlib import Path

from scalar_convergence_gate import available_times, evaluate_contract

def read(p):
    p = Path(p)
    return p.read_text() if p.exists() else ""

def field_average_block(control_dict):
    if "fieldAverage" not in control_dict:
        return ""
    start = control_dict.find("fieldAverage")
    return control_dict[start:start + 1200]

def requests_tprime(control_dict):
    fa = field_average_block(control_dict)
    return bool(re.search(r"\bT\s*\{[^{}]*\bprime2Mean\s+on\b", fa, re.S))

def wall_temperature_series(case):
    rows = []
    for path in sorted((case / "postProcessing" / "fwT").glob("*/surfaceFieldValue.dat")):
        for line in path.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                rows.append((float(parts[0]), float(parts[-1])))
            except ValueError:
                continue
    # De-duplicate repeated times if a function object got restarted.
    dedup = {}
    for t, val in rows:
        dedup[t] = val
    return sorted(dedup.items())

def half_window_drift_k(series):
    if len(series) < 6:
        return None
    vals = [v for _, v in series]
    mid = len(vals) // 2
    if mid == 0 or mid == len(vals):
        return None
    first = vals[:mid]
    second = vals[mid:]
    return abs(sum(second) / len(second) - sum(first) / len(first))

def geometry_avg_threshold(case):
    meta = case / "geometry_c2.json"
    if not meta.exists():
        return None
    try:
        data = json.loads(meta.read_text())
        params = data.get("parameters", data.get("params", data))
        radius = float(params["radius"])
        arc_deg = float(params["arc_deg"])
        u_in = float(params["u_in"])
    except Exception:
        return None
    if u_in <= 0:
        return None
    ft = math.radians(arc_deg) * radius / u_in
    return 3.0 * ft, ft, u_in

def geometry_label_u(case):
    meta = case / "geometry_c2.json"
    if not meta.exists():
        return None
    try:
        data = json.loads(meta.read_text())
        params = data.get("parameters", data.get("params", data))
        return float(params["u_in"])
    except Exception:
        return None

def boundary_data_bulk(case, patch="inlet"):
    path = case / "constant" / "boundaryData" / patch / "0" / "U"
    if not path.exists():
        return None
    vecs = [
        [float(v) for v in n.split()]
        for n in re.findall(r"\(([^()]*)\)", path.read_text())
    ]
    if not vecs:
        return None
    mags = [math.sqrt(x * x + y * y + z * z) for x, y, z in vecs]
    return sum(mags) / len(mags)

def main(case):
    case = Path(case)
    cd = read(case/"system"/"controlDict")
    fs = read(case/"system"/"fvSolution")
    dp = read(case/"system"/"decomposeParDict")
    cm = read(case/"log.checkMesh")
    results = []  # (gotcha, status, detail)
    def add(g, ok, detail, warn=False):
        results.append((g, "PASS" if ok else ("WARN" if warn else "FAIL"), detail))

    # #1 physical-second times (heuristic: deltaT present and < 0.01 s, endTime a plausible physical time)
    m = re.search(r"\bdeltaT\s+([0-9.eE+-]+)", cd)
    dt = float(m.group(1)) if m else None
    add("1 physical-second times", dt is not None and dt < 0.05,
        f"deltaT={dt}" if dt else "deltaT not found", warn=(dt is not None and dt>=0.05))

    # #2 pRefCell/pRefValue (needed for fully-periodic/closed incompressible; N/A if a fixedValue p patch pins it)
    has_pref = "pRefCell" in fs and "pRefValue" in fs
    pfield = read(case/"0"/"p")
    p_pinned = bool(re.search(r"type\s+fixedValue", pfield))  # open inlet/outlet case fixes p at a boundary
    if p_pinned and not has_pref:
        add("2 pRefCell/pRefValue", True, "N/A: open case, fixedValue p patch pins the level")
    else:
        add("2 pRefCell/pRefValue", has_pref, "present in fvSolution" if has_pref else "MISSING (periodic case will fail)")

    # #3 UFinal solver entry (explicit or via regex)
    has_ufinal = bool(re.search(r'UFinal|"[^"]*U[^"]*Final[^"]*"', fs))
    add("3 UFinal solver entry", has_ufinal, "found UFinal/regex" if has_ufinal else "MISSING UFinal solver entry")

    # #4 fieldAverage writeControl writeTime + no per-step bloat
    if "fieldAverage" in cd:
        fa = field_average_block(cd)
        wt = "writeControl" in fa and "writeTime" in fa[:fa.find("}")+400]
        add("4 fieldAverage write-control", wt, "writeControl writeTime" if wt else "fieldAverage not writeTime -> bloat risk")
        ts = re.search(r"\btimeStart\s+([0-9.eE+-]+)", fa)
        time_start = float(ts.group(1)) if ts else None
        threshold = geometry_avg_threshold(case)
        if time_start is None:
            add("4b fieldAverage start hygiene", False, "fieldAverage present but timeStart not found")
        elif threshold:
            min_start, ft, u_meta = threshold
            add(
                "4b fieldAverage start hygiene",
                time_start >= min_start,
                f"timeStart={time_start:.6g}; required >= ramp+1FT ~= 3FT = {min_start:.6g} s "
                f"(FT={ft:.6g} s from geometry metadata u_in={u_meta:.6g})",
            )
        else:
            add(
                "4b fieldAverage start hygiene",
                True,
                f"timeStart={time_start:.6g}; geometry metadata unavailable, cannot verify ramp+1FT",
                warn=True,
            )
    else:
        add("4 fieldAverage write-control", True, "no fieldAverage (n/a)")
        add("4b fieldAverage start hygiene", True, "no fieldAverage (n/a)")

    if requests_tprime(cd):
        series = wall_temperature_series(case)
        drift = half_window_drift_k(series)
        if drift is None:
            add(
                "4d variance-stat replicate-window gate",
                False,
                "TPrime2Mean requested; no usable fwT half-window series yet. Before promoting variance/flicker, "
                "use >=2 independent developed averaging windows and report the spread; one window is not enough.",
                warn=True,
            )
        else:
            add(
                "4d variance-stat replicate-window gate",
                False,
                f"TPrime2Mean requested; fwT half-window areaAverage(T) drift={drift:.3g} K "
                "(diagnostic only, not a pass/fail threshold). Promote variance/flicker only with >=2 "
                "independent developed windows and reported scatter; warning does not block mean-field setup.",
                warn=True,
            )
    else:
        add("4d variance-stat replicate-window gate", True, "TPrime2Mean not requested (n/a)")

    # #4e exact source balance and bulk-field flatness. A source-driven scalar
    # must carry an explicit balance contract, and completed output must pass it
    # before any profile/statistic can be promoted.
    has_scalar_source = "scalarSemiImplicitSource" in cd
    contract_path = case / "constant" / "sourceDrivenScalarGate.json"
    if has_scalar_source:
        if not contract_path.exists():
            add(
                "4e source-driven scalar equilibrium",
                False,
                "scalarSemiImplicitSource present but sourceDrivenScalarGate.json is missing; "
                "exact source/sink balance is a hard production requirement",
            )
        else:
            try:
                contract = json.loads(contract_path.read_text())
                first = contract["scalars"][0]
                times = available_times(case, first["field"])
                start = float(contract["flatness_time_start"])
                minimum = int(contract.get("minimum_flatness_snapshots", 3))
                developed = [row for row in times if row[0] + 1e-12 >= start]
                if len(developed) < minimum:
                    add(
                        "4e source-driven scalar equilibrium",
                        False,
                        f"balance contract configured; runtime gate pending "
                        f"({len(developed)}/{minimum} snapshots at t>={start:g})",
                        warn=True,
                    )
                else:
                    gate = evaluate_contract(case, contract_path)
                    worst_balance = max(
                        row["balance_relative_error"] for row in gate["scalars"]
                    )
                    worst_flatness = max(
                        row["flatness_relative_range"] for row in gate["scalars"]
                    )
                    add(
                        "4e source-driven scalar equilibrium",
                        gate["statistics_allowed"],
                        f"{gate['status']}: worst source-balance error={100*worst_balance:.2f}%, "
                        f"bulk range={100*worst_flatness:.3f}%",
                    )
            except Exception as exc:
                add(
                    "4e source-driven scalar equilibrium",
                    False,
                    f"could not evaluate balance contract: {exc}",
                )
    else:
        add("4e source-driven scalar equilibrium", True, "no active scalar source (n/a)")

    # Momentum/bulk flatness check for meanVelocityForce channels. The forcing
    # log supplies a reference-free delivered-Ubar history.
    if "meanVelocityForce" in read(case / "constant" / "fvOptions"):
        log_text = read(case / "log.pimpleFoam")
        ubar = [
            float(value)
            for value in re.findall(r"uncorrected Ubar\s*=\s*([0-9.eE+-]+)", log_text)
        ]
        if len(ubar) >= 20:
            tail = ubar[-min(200, len(ubar)):]
            relative_range = (max(tail) - min(tail)) / max(abs(sum(tail) / len(tail)), 1e-30)
            add(
                "4f momentum bulk flatness",
                relative_range <= 0.01,
                f"last {len(tail)} delivered-Ubar samples span {100*relative_range:.4f}% "
                "(hard limit 1%)",
            )
        else:
            add(
                "4f momentum bulk flatness",
                False,
                "meanVelocityForce present; delivered-Ubar history unavailable/pending",
                warn=True,
            )
    else:
        add("4f momentum bulk flatness", True, "meanVelocityForce not used (n/a)")

    label_u = geometry_label_u(case)
    delivered_u = boundary_data_bulk(case)
    if label_u is not None and delivered_u is not None and label_u > 0:
        mismatch = abs(delivered_u - label_u) / label_u
        add(
            "4c operating-point label vs delivered bulk",
            mismatch <= 0.10,
            f"geometry label u_in={label_u:.6g}, boundaryData mean|U|={delivered_u:.6g}, mismatch={mismatch*100:.1f}%",
            warn=True,
        )
    elif delivered_u is not None:
        add(
            "4c operating-point label vs delivered bulk",
            True,
            f"boundaryData mean|U|={delivered_u:.6g}; geometry label unavailable",
            warn=True,
        )
    else:
        add("4c operating-point label vs delivered bulk", True, "no mapped boundaryData inlet (n/a)")

    # #5 numberOfSubdomains in decomposeParDict
    m = re.search(r"numberOfSubdomains\s+(\d+)", dp)
    nsub = int(m.group(1)) if m else None
    add("5 numberOfSubdomains in-dict", nsub is not None, f"numberOfSubdomains={nsub}" if nsub else "MISSING (the -flag is unsupported in 2512)")

    # #6 decomposition method. scotch is the RECOMMENDED default for curved geometry tiers (OPT-3:
    # `simple` load-imbalances 3.92x on a non-box mesh); it is valid as long as libscotchDecomp.so is
    # present (verified locally + standard on Grace). simple/hierarchical are valid for box channels.
    method = re.search(r"\bmethod\s+(\w+)", dp)
    meth = method.group(1) if method else None
    scotch_ok = bool(glob.glob(str(Path(os.environ.get("FOAM_LIBBIN", "/usr/lib/openfoam/openfoam2512/platforms/linux64GccDPInt32Opt/lib")) / "libscotchDecomp.so")))
    if meth == "scotch":
        add("6 decomposition method", scotch_ok, "method=scotch (recommended for geometry)" if scotch_ok
            else "method=scotch but libscotchDecomp.so NOT found -> falls back/fails", warn=not scotch_ok)
    else:
        add("6 decomposition method", meth in ("simple","hierarchical"), f"method={meth}" if meth else "no method")

    # cells/core band (25k-50k) -- needs cell count from checkMesh
    cm_cells = re.search(r"cells:\s+(\d+)", cm)
    ncells = int(cm_cells.group(1)) if cm_cells else None
    if ncells and nsub:
        cpc = ncells/nsub
        add("6b cells/core 25k-50k", 20000 <= cpc <= 60000, f"{ncells} cells / {nsub} = {cpc:.0f} cells/core",
            warn=not (20000 <= cpc <= 60000))
    else:
        add("6b cells/core 25k-50k", True, "cell count or nsub unknown (run checkMesh)", warn=True)

    # #7 mpirun (runtime-only reminder)
    add("7 mpirun root/slots", True, "RUNTIME: --allow-run-as-root + OMPI_ALLOW_RUN_AS_ROOT[_CONFIRM]=1; "
        "OpenMPI counts PHYSICAL cores (use --use-hwthread-cpus); on Grace use srun", warn=True)

    # checkMesh OK?
    add("mesh", "Mesh OK" in cm, "checkMesh: Mesh OK" if "Mesh OK" in cm else "checkMesh not OK / not run", warn=("Mesh OK" not in cm))

    print(f"=== HPC pre-flight: {case} ===")
    nfail = 0
    for g, st, d in results:
        mark = {"PASS":"[ OK ]","WARN":"[WARN]","FAIL":"[FAIL]"}[st]
        print(f"  {mark} gotcha {g}: {d}")
        if st == "FAIL": nfail += 1
    print(f"--- {nfail} FAIL ---" if nfail else "--- all gotchas clear ---")
    return 1 if nfail else 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
