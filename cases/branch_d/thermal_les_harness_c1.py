#!/usr/bin/env python3
"""C1 thermal LES postprocessor.

Reads OpenFOAM ASCII fields from a channel LES case and writes a wall-normal
thermal profile.  The harness is intentionally conservative: without local
Kawamura profile data it reports "profile generated" instead of pass/fail.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")


def read_text(path: Path) -> str:
    return path.read_text(errors="replace")


def parse_nonuniform_block(text: str, field_type: str) -> str:
    pattern = rf"internalField\s+nonuniform\s+List<{field_type}>\s*\d+\s*\((.*?)\)\s*;"
    match = re.search(pattern, text, re.S)
    if not match:
        raise ValueError(f"no nonuniform List<{field_type}> internalField")
    return match.group(1)


def read_scalar(path: Path) -> np.ndarray:
    block = parse_nonuniform_block(read_text(path), "scalar")
    return np.array([float(x) for x in block.split()], dtype=float)


def read_vector(path: Path) -> np.ndarray:
    block = parse_nonuniform_block(read_text(path), "vector")
    rows = re.findall(r"\(([^()]*)\)", block)
    return np.array([[float(v) for v in row.split()] for row in rows], dtype=float)


def numeric_times(case: Path) -> list[str]:
    out: list[tuple[float, str]] = []
    for entry in case.iterdir():
        if not entry.is_dir():
            continue
        try:
            out.append((float(entry.name), entry.name))
        except ValueError:
            pass
    return [name for _, name in sorted(out)]


def accumulated_time(case: Path, time_name: str, field: str) -> float:
    path = case / time_name / "uniform/functionObjects/functionObjectProperties"
    if not path.exists():
        raise FileNotFoundError(path)
    pattern = re.compile(
        rf"(?ms)^\s*{re.escape(field)}\s*\{{.*?^\s*totalTime\s+([\d.eE+-]+);"
    )
    match = pattern.search(read_text(path))
    if not match:
        raise ValueError(f"no fieldAverage totalTime for {field} at {time_name}")
    return float(match.group(1))


def accumulated_window_mean(
    case: Path,
    opening_time: str,
    endpoint_time: str,
    mean_field: str,
    base_field: str,
    reader,
) -> tuple[np.ndarray, list[float]]:
    early_total = accumulated_time(case, opening_time, base_field)
    late_total = accumulated_time(case, endpoint_time, base_field)
    if late_total <= early_total:
        raise ValueError(
            f"non-increasing accumulation time for {base_field}: "
            f"{early_total} -> {late_total}"
        )
    early = reader(case / opening_time / mean_field)
    late = reader(case / endpoint_time / mean_field)
    return (
        (late * late_total - early * early_total) / (late_total - early_total),
        [early_total, late_total],
    )


def latest_time_with(case: Path, required: list[str]) -> str:
    for time_name in reversed(numeric_times(case)):
        tdir = case / time_name
        if all((tdir / field).exists() for field in required):
            return time_name
    raise FileNotFoundError(f"no time directory has required fields: {required}")


def matching_static_scalar(
    case: Path, preferred_time: str, field: str, expected_length: int
) -> tuple[np.ndarray, Path]:
    """Load a static field without forcing the analysis back to its write time."""
    candidates = [preferred_time]
    candidates.extend(reversed(numeric_times(case)))
    seen: set[str] = set()
    for time_name in candidates:
        if time_name in seen:
            continue
        seen.add(time_name)
        path = case / time_name / field
        if not path.exists():
            continue
        values = read_scalar(path)
        if len(values) == expected_length:
            return values, path
    raise FileNotFoundError(
        f"no root {field} field with {expected_length} values under {case}"
    )


def parse_last_pressure_gradient(log: Path) -> float | None:
    if not log.exists():
        return None
    vals = re.findall(r"pressure gradient =\s*([-+0-9.eE]+)", read_text(log))
    return float(vals[-1]) if vals else None


def pressure_gradient_window(log: Path, start: float, end: float | None = None) -> np.ndarray:
    """Return one final pressure-gradient value per solver time at/after start."""
    if not log.exists():
        return np.asarray([], dtype=float)
    current_time: float | None = None
    by_time: dict[float, float] = {}
    for line in read_text(log).splitlines():
        time_match = re.match(r"\s*Time\s*=\s*([-+0-9.eE]+)\s*$", line)
        if time_match:
            current_time = float(time_match.group(1))
            continue
        gradient_match = re.search(r"pressure gradient =\s*([-+0-9.eE]+)", line)
        if (
            gradient_match
            and current_time is not None
            and current_time >= start
            and (end is None or current_time <= end)
        ):
            by_time[current_time] = float(gradient_match.group(1))
    return np.asarray([by_time[key] for key in sorted(by_time)], dtype=float)


def grouped_profile(cy: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ylevels = np.unique(np.round(cy, 10))
    y = []
    prof = []
    for yl in ylevels:
        mask = np.isclose(cy, yl)
        y.append(cy[mask].mean())
        prof.append(values[mask].mean(axis=0) if values.ndim == 2 else values[mask].mean())
    return np.array(y), np.array(prof)


def y_groups(cy: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rounded = np.round(cy, 10)
    ylevels, inverse = np.unique(rounded, return_inverse=True)
    counts = np.bincount(inverse).astype(float)
    y = np.bincount(inverse, weights=cy) / counts
    return y, inverse, counts


def snapshot_times_with(
    case: Path,
    required: list[str],
    start: float,
    end: float | None,
    max_snapshots: int,
) -> list[str]:
    valid: list[tuple[float, str]] = []
    for time_name in numeric_times(case):
        try:
            tval = float(time_name)
        except ValueError:
            continue
        if tval < start:
            continue
        if end is not None and tval > end:
            continue
        tdir = case / time_name
        if all((tdir / field).exists() for field in required):
            valid.append((tval, time_name))
    if max_snapshots > 0 and len(valid) > max_snapshots:
        valid = valid[-max_snapshots:]
    return [name for _, name in valid]


def compute_snapshot_covariance(
    case: Path,
    cy: np.ndarray,
    scalar_field: str,
    start: float,
    end: float | None,
    max_snapshots: int,
    nu: float,
    alpha: float,
) -> tuple[list[str], dict[str, np.ndarray]]:
    times = snapshot_times_with(
        case,
        ["U", scalar_field],
        start=start,
        end=end,
        max_snapshots=max_snapshots,
    )
    if not times:
        raise FileNotFoundError(f"no U/{scalar_field} snapshots available for covariance")

    y, inverse, _ = y_groups(cy)
    nlev = len(y)
    accum = {name: np.zeros(nlev) for name in ("n", "u", "v", "t", "uv", "vT", "t2")}

    for time_name in times:
        tdir = case / time_name
        u = read_vector(tdir / "U")
        temp = read_scalar(tdir / scalar_field)
        if len(temp) != len(cy) or len(u) != len(cy):
            raise ValueError(f"{time_name}: field sizes do not match Cy")
        ux = u[:, 0]
        uy = u[:, 1]
        accum["n"] += np.bincount(inverse, minlength=nlev)
        accum["u"] += np.bincount(inverse, weights=ux, minlength=nlev)
        accum["v"] += np.bincount(inverse, weights=uy, minlength=nlev)
        accum["t"] += np.bincount(inverse, weights=temp, minlength=nlev)
        accum["uv"] += np.bincount(inverse, weights=ux * uy, minlength=nlev)
        accum["vT"] += np.bincount(inverse, weights=uy * temp, minlength=nlev)
        accum["t2"] += np.bincount(inverse, weights=temp * temp, minlength=nlev)

    n = accum["n"]
    u_mean = accum["u"] / n
    v_mean = accum["v"] / n
    t_mean = accum["t"] / n
    uv_cov = accum["uv"] / n - u_mean * v_mean
    vT_cov = accum["vT"] / n - v_mean * t_mean
    t_var = np.maximum(accum["t2"] / n - t_mean * t_mean, 0.0)

    order = np.argsort(y)
    y = y[order]
    u_mean = u_mean[order]
    t_mean = t_mean[order]
    uv_cov = uv_cov[order]
    vT_cov = vT_cov[order]
    t_var = t_var[order]

    dudy = finite_gradient(y, u_mean)
    dtdy = finite_gradient(y, t_mean)
    with np.errstate(divide="ignore", invalid="ignore"):
        nu_t = -uv_cov / dudy
        alpha_t = -vT_cov / dtdy
        pr_t = nu_t / alpha_t
        nu_t_plus = nu_t / nu
        alpha_t_plus = alpha_t / alpha

    finite_positive = (nu_t > 0) & (alpha_t > 0) & np.isfinite(pr_t)
    pr_t = np.where(finite_positive, pr_t, np.nan)

    stats = {
        "y": y,
        "U_mean_snapshot": u_mean,
        "T_mean_snapshot": t_mean,
        "uv_resolved": uv_cov,
        "vT_resolved": vT_cov,
        "T_rms_snapshot": np.sqrt(t_var),
        "dUdy_snapshot": dudy,
        "dTdy_snapshot": dtdy,
        "nu_t_resolved": nu_t,
        "alpha_t_resolved": alpha_t,
        "nu_t_over_nu": nu_t_plus,
        "alpha_t_over_alpha": alpha_t_plus,
        "Pr_t_resolved": pr_t,
    }
    return times, stats


def finite_gradient(y: np.ndarray, f: np.ndarray) -> np.ndarray:
    return np.gradient(f, y, edge_order=1)


def write_profile_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_profile(path: Path, yp: np.ndarray, tplus: np.ndarray, trms_plus: np.ndarray | None) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional plotting dependency
        print(f"plot skipped: {exc}")
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.semilogx(yp, tplus, "o-", ms=4, color="#3B6EA8", label="T_plus profile")
    if trms_plus is not None:
        ax.semilogx(yp, trms_plus, "s-", ms=3, color="#B85C38", label="T_rms_plus")
    ax.set_xlabel("y_plus")
    ax.set_ylabel("thermal wall units")
    ax.set_title("C1 thermal LES profile scaffold")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def plot_prt(path: Path, yp: np.ndarray, pr_t: np.ndarray) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional plotting dependency
        print(f"Pr_t plot skipped: {exc}")
        return

    good = np.isfinite(pr_t)
    if not np.any(good):
        print("Pr_t plot skipped: no finite positive eddy-diffusivity points")
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.semilogx(yp[good], pr_t[good], "o-", ms=4, color="#4C7A3F", label="resolved Pr_t")
    ax.axhspan(0.7, 1.0, color="grey", alpha=0.16, label="typical turbulent Pr_t band")
    ax.set_xlabel("y_plus")
    ax.set_ylabel("resolved turbulent Pr_t")
    ax.set_title("C1 snapshot covariance scaffold")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", nargs="?", default="les_channel_c1_thermal_smoke")
    parser.add_argument("time", nargs="?")
    parser.add_argument("--nu", type=float, default=3.09278e-6)
    parser.add_argument("--pr", type=float, default=14.4)
    parser.add_argument("--delta", type=float, default=0.01)
    parser.add_argument("--field", default="T", help="thermal scalar field to read, e.g. T or H")
    parser.add_argument("--dpdx", type=float, help="kinematic pressure gradient magnitude [m/s^2]")
    parser.add_argument(
        "--solver-log",
        default="log.pimpleFoam",
        help="solver log path, absolute or relative to the case directory",
    )
    parser.add_argument(
        "--dpdx-window-start",
        type=float,
        help="average one final pressure-gradient value per solver time from this time onward",
    )
    parser.add_argument("--dpdx-window-end", type=float)
    parser.add_argument(
        "--mean-window-start",
        type=float,
        help=(
            "derive UMean/scalarMean over the late window using cumulative "
            "fieldAverage writes; selects the first write at/after this time"
        ),
    )
    parser.add_argument(
        "--wall-value",
        type=float,
        help="known lower-wall scalar value; enables a wall-anchored gradient fit",
    )
    parser.add_argument(
        "--wall-fit-levels",
        type=int,
        default=3,
        help="near-wall cell-center levels used with --wall-value (default: 3)",
    )
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--covariance", action="store_true", help="compute resolved uv, vT, and Pr_t from saved U/T snapshots")
    parser.add_argument("--covariance-start", type=float, default=0.0, help="first snapshot time for covariance")
    parser.add_argument("--covariance-end", type=float, help="last snapshot time for covariance")
    parser.add_argument("--max-snapshots", type=int, default=80, help="limit covariance to the latest N snapshots; <=0 uses all")
    args = parser.parse_args()

    case = Path(args.case)
    solver_log = Path(args.solver_log)
    if not solver_log.is_absolute():
        solver_log = case / solver_log
    alpha = args.nu / args.pr
    if args.time:
        time_name = args.time
    else:
        scalar_mean_field = f"{args.field}Mean"
        try:
            time_name = latest_time_with(case, [scalar_mean_field])
        except FileNotFoundError:
            time_name = latest_time_with(case, [args.field])
    tdir = case / time_name

    scalar_mean_field = f"{args.field}Mean"
    scalar_prime_field = f"{args.field}Prime2Mean"
    u_field = "UMean" if (tdir / "UMean").exists() else "U"
    t_field = scalar_mean_field if (tdir / scalar_mean_field).exists() else args.field
    mean_window = None
    trms_values = None
    if args.mean_window_start is not None:
        if u_field != "UMean" or t_field != scalar_mean_field:
            raise ValueError("--mean-window-start requires UMean and scalarMean fields")
        eligible = [
            name
            for name in numeric_times(case)
            if float(name) + 1e-12 >= args.mean_window_start
            and float(name) < float(time_name) - 1e-12
            and (case / name / "UMean").exists()
            and (case / name / scalar_mean_field).exists()
        ]
        if not eligible:
            raise ValueError(
                f"no accumulator write opens mean window at/after {args.mean_window_start}"
            )
        opening_time = eligible[0]
        u, u_totals = accumulated_window_mean(
            case, opening_time, time_name, "UMean", "U", read_vector
        )
        temp, scalar_totals = accumulated_window_mean(
            case, opening_time, time_name, scalar_mean_field, args.field, read_scalar
        )
        if max(abs(a - b) for a, b in zip(u_totals, scalar_totals)) > 1e-10:
            raise ValueError(
                f"U/scalar accumulation times differ: {u_totals} vs {scalar_totals}"
            )
        mean_window = {
            "requested_start": args.mean_window_start,
            "opening_time": opening_time,
            "endpoint_time": time_name,
            "accumulation_times": scalar_totals,
        }
        u_field = f"UMean window [{opening_time}, {time_name}]"
        t_field = f"{scalar_mean_field} window [{opening_time}, {time_name}]"
        if (case / opening_time / scalar_prime_field).exists() and (
            tdir / scalar_prime_field
        ).exists():
            early_mean = read_scalar(case / opening_time / scalar_mean_field)
            late_mean = read_scalar(tdir / scalar_mean_field)
            early_var = read_scalar(case / opening_time / scalar_prime_field)
            late_var = read_scalar(tdir / scalar_prime_field)
            early_total, late_total = scalar_totals
            raw_window = (
                (late_var + late_mean**2) * late_total
                - (early_var + early_mean**2) * early_total
            ) / (late_total - early_total)
            trms_values = np.sqrt(np.clip(raw_window - temp**2, 0, None))
    else:
        u = read_vector(tdir / u_field)
        temp = read_scalar(tdir / t_field)
    if len(u) != len(temp):
        raise ValueError(
            f"{time_name}: {u_field} has {len(u)} cells but {t_field} has {len(temp)}"
        )
    cy, cy_path = matching_static_scalar(case, time_name, "Cy", len(temp))

    y, uprof = grouped_profile(cy, u[:, 0])
    _, tprof = grouped_profile(cy, temp)
    order = np.argsort(y)
    y = y[order]
    uprof = uprof[order]
    tprof = tprof[order]

    dpdx_samples = np.asarray([], dtype=float)
    if args.dpdx is not None:
        dpdx = args.dpdx
        dpdx_method = "explicit --dpdx"
    elif args.dpdx_window_start is not None:
        dpdx_samples = pressure_gradient_window(
            solver_log, args.dpdx_window_start, args.dpdx_window_end
        )
        if not len(dpdx_samples):
            raise ValueError(
                f"no pressure-gradient samples at/after {args.dpdx_window_start}"
            )
        dpdx = float(dpdx_samples.mean())
        end_tag = (
            f" to t<={args.dpdx_window_end:g}"
            if args.dpdx_window_end is not None
            else ""
        )
        dpdx_method = (
            f"window mean from {solver_log} at t>={args.dpdx_window_start:g}{end_tag}"
        )
    else:
        dpdx = parse_last_pressure_gradient(solver_log)
        dpdx_method = f"last sample from {solver_log}"
    if dpdx is None:
        dpdx = 0.46917746
        dpdx_method = "fallback default O pressure gradient"
        print(f"warning: using default O pressure gradient {dpdx:g}")
    utau = math.sqrt(abs(dpdx) * args.delta)
    retau = utau * args.delta / args.nu

    dtdy = finite_gradient(y, tprof)
    if args.wall_value is not None:
        fit_levels = min(args.wall_fit_levels, len(y))
        if fit_levels < 1:
            raise ValueError("--wall-fit-levels must be at least 1")
        wall_grad = float(
            np.polyfit(
                np.r_[0.0, y[:fit_levels]],
                np.r_[args.wall_value, tprof[:fit_levels]],
                1,
            )[0]
        )
        twall = args.wall_value
        wall_gradient_method = f"fixed-wall {fit_levels}-level linear fit"
    else:
        wall_grad = dtdy[0]
        twall = tprof[0] - wall_grad * y[0]
        wall_gradient_method = "first-cell finite difference with extrapolated wall"
    ttau = abs(alpha * wall_grad) / utau if abs(wall_grad) > 0 else float("nan")

    ywall = np.minimum(y, 2 * args.delta - y)
    yp = ywall * utau / args.nu
    tplus = np.abs(twall - tprof) / ttau if np.isfinite(ttau) and ttau > 0 else np.full_like(tprof, np.nan)

    trms_plus = None
    if trms_values is None and (tdir / scalar_prime_field).exists():
        trms_values = np.sqrt(np.clip(read_scalar(tdir / scalar_prime_field), 0, None))
    if trms_values is not None and np.isfinite(ttau) and ttau > 0:
        _, trms_prof = grouped_profile(cy, trms_values)
        trms_plus = trms_prof[order] / ttau

    half = y <= args.delta
    half_order = np.argsort(yp[half])
    yp_h = yp[half][half_order]
    tplus_h = tplus[half][half_order]
    trms_h = trms_plus[half][half_order] if trms_plus is not None else None
    y_h = ywall[half][half_order]
    u_h = uprof[half][half_order]
    t_h = tprof[half][half_order]

    outdir = Path(args.outdir) if args.outdir else case.parent / "figs"
    stem = case.name
    rows: list[dict[str, float]] = []
    covariance_stats = None
    covariance_times: list[str] = []
    if args.covariance:
        covariance_times, covariance_stats = compute_snapshot_covariance(
            case=case,
            cy=cy,
            scalar_field=args.field,
            start=args.covariance_start,
            end=args.covariance_end,
            max_snapshots=args.max_snapshots,
            nu=args.nu,
            alpha=alpha,
        )
        print(
            f"covariance snapshots: {len(covariance_times)} "
            f"from {covariance_times[0]} to {covariance_times[-1]}"
        )

    for i in range(len(yp_h)):
        row = {
            "y": float(y_h[i]),
            "y_plus": float(yp_h[i]),
            "U_mean": float(u_h[i]),
            "T_mean": float(t_h[i]),
            "T_plus": float(tplus_h[i]),
        }
        if trms_h is not None:
            row["T_rms_plus"] = float(trms_h[i])
        if covariance_stats is not None:
            cov_y = covariance_stats["y"]
            j = int(np.argmin(np.abs(cov_y - y_h[i])))
            row.update(
                {
                    "uv_resolved": float(covariance_stats["uv_resolved"][j]),
                    "vT_resolved": float(covariance_stats["vT_resolved"][j]),
                    "dUdy_snapshot": float(covariance_stats["dUdy_snapshot"][j]),
                    "dTdy_snapshot": float(covariance_stats["dTdy_snapshot"][j]),
                    "nu_t_resolved": float(covariance_stats["nu_t_resolved"][j]),
                    "alpha_t_resolved": float(covariance_stats["alpha_t_resolved"][j]),
                    "nu_t_over_nu": float(covariance_stats["nu_t_over_nu"][j]),
                    "alpha_t_over_alpha": float(covariance_stats["alpha_t_over_alpha"][j]),
                    "Pr_t_resolved": float(covariance_stats["Pr_t_resolved"][j]),
                    "T_rms_snapshot": float(covariance_stats["T_rms_snapshot"][j]),
                }
            )
        rows.append(row)

    csv_path = outdir / f"c1_thermal_profile_{stem}.csv"
    write_profile_csv(csv_path, rows)
    print(f"=== C1 thermal harness: {case}/{time_name} ===")
    print(f"fields: U={u_field} scalar={t_field} base_scalar={args.field}")
    print(f"coordinate provenance: {cy_path}")
    print(f"Pr={args.pr:g} alpha={alpha:.6g} u_tau={utau:.6g} Re_tau={retau:.1f}")
    print(
        f"pressure-gradient provenance: {dpdx_method}; value={dpdx:.9g}; "
        f"samples={len(dpdx_samples)}"
    )
    print(f"near-wall dT/dy={wall_grad:.6g} T_tau={ttau:.6g}")
    print(f"wall normalization: {wall_gradient_method}; wall value={twall:.9g}")
    print(f"wrote {csv_path}")
    print("Kawamura comparison: pending local DNS profile/digitization; no pass/fail assigned.")

    tplus_png = outdir / f"c1_thermal_Tplus_{stem}.png"
    prt_png = outdir / f"c1_thermal_Prt_{stem}.png"
    plot_profile(tplus_png, yp_h, tplus_h, trms_h)
    prt_written = False
    if covariance_stats is not None:
        cov_y = covariance_stats["y"]
        ywall_cov = np.minimum(cov_y, 2 * args.delta - cov_y)
        yp_cov = ywall_cov * utau / args.nu
        half_cov = cov_y <= args.delta
        order_cov = np.argsort(yp_cov[half_cov])
        plot_prt(
            prt_png,
            yp_cov[half_cov][order_cov],
            covariance_stats["Pr_t_resolved"][half_cov][order_cov],
        )
        prt_written = prt_png.exists()

    summary = {
        "case": str(case),
        "case_name": case.name,
        "time": time_name,
        "fields": {"U": u_field, "thermal_scalar": t_field, "base_thermal_scalar": args.field},
        "coordinate_field": {"name": "Cy", "path": str(cy_path), "length": len(cy)},
        "nu": args.nu,
        "Pr": args.pr,
        "alpha": alpha,
        "delta": args.delta,
        "dpdx_kinematic": dpdx,
        "dpdx_method": dpdx_method,
        "solver_log": str(solver_log),
        "mean_window": mean_window,
        "dpdx_sample_count": len(dpdx_samples),
        "dpdx_sample_std": float(dpdx_samples.std()) if len(dpdx_samples) else None,
        "u_tau": utau,
        "Re_tau": retau,
        "wall_dTdy": wall_grad,
        "T_tau": ttau,
        "T_wall_linear_extrapolated": twall,
        "wall_gradient_method": wall_gradient_method,
        "wall_value_input": args.wall_value,
        "wall_fit_levels": args.wall_fit_levels if args.wall_value is not None else None,
        "profile_rows": len(rows),
        "y_plus_min": float(np.nanmin(yp_h)) if len(yp_h) else None,
        "y_plus_max": float(np.nanmax(yp_h)) if len(yp_h) else None,
        "outputs": {
            "profile_csv": str(csv_path),
            "Tplus_png": str(tplus_png),
            "Prt_png": str(prt_png) if prt_written else None,
        },
        "covariance": {
            "enabled": covariance_stats is not None,
            "snapshot_count": len(covariance_times),
            "first_snapshot": covariance_times[0] if covariance_times else None,
            "last_snapshot": covariance_times[-1] if covariance_times else None,
            "requested_start": args.covariance_start,
            "requested_end": args.covariance_end,
            "max_snapshots": args.max_snapshots,
        },
        "certification_status": "machinery/scaffold only; no DNS pass/fail assigned",
    }
    summary_path = outdir / f"c1_thermal_summary_{stem}.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
