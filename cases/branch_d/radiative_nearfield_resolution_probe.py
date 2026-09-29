#!/usr/bin/env python3
"""Diagnostic probe for the withdrawn v2 line-of-sight radiative spike.

This is not a replacement radiation model and must not be used to promote a v2
hotspot number. It asks a narrower numerical question: how sensitive is the
withdrawn radiative peak/mean to enforcing a minimum wall-resolvable length?
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


LIVE_ROOT = Path(__file__).resolve().parents[2]
IN_LIVE_TREE = (LIVE_ROOT / "cases" / "branch_d").is_dir()
DEFAULT_IN = (
    LIVE_ROOT / "results/advisor_packet/99_plasma_qmap_v2.json"
    if IN_LIVE_TREE else Path("figs/plasma_qmap_v2.json")
)
DEFAULT_OUT = (
    LIVE_ROOT / "cases/branch_d/figs" if IN_LIVE_TREE else Path("figs")
)
DEFAULT_ADVISOR_OUT = LIVE_ROOT / "results/advisor_packet" if IN_LIVE_TREE else None


def periodic_distance_matrix(s: np.ndarray, period: float) -> np.ndarray:
    dx = np.abs(s[:, None] - s[None, :])
    return np.minimum(dx, period - dx)


def gaussian_periodic_smooth(y: np.ndarray, s: np.ndarray, period: float, fwhm: float) -> np.ndarray:
    if fwhm <= 0:
        out = y.copy()
    else:
        sigma = fwhm / 2.354820045
        dist = periodic_distance_matrix(s, period)
        w = np.exp(-0.5 * (dist / sigma) ** 2)
        w /= w.sum(axis=1, keepdims=True)
        out = w @ y
    return out / out.mean()


def contiguous_fwhm_width(s: np.ndarray, y: np.ndarray, period: float) -> float:
    y = y / y.mean()
    imax = int(np.argmax(y))
    half = 0.5 * (1.0 + float(y[imax]))
    order = np.r_[np.arange(imax, len(y)), np.arange(0, imax)]
    ss = np.r_[s[imax:], s[:imax] + period]
    yy = y[order]
    above = yy >= half
    if not above.any():
        return 0.0
    idx = np.where(above)[0]
    return float(ss[idx[-1]] - ss[idx[0]] + (s[1] - s[0]))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=DEFAULT_IN)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--advisor-out", type=Path, default=DEFAULT_ADVISOR_OUT)
    args = ap.parse_args()

    data = json.loads(args.input.read_text())
    s = np.asarray(data["s_m"], dtype=float)
    ds = float(np.median(np.diff(s)))
    period = float(s[-1] + ds)
    q_rad = np.asarray(data["q_rad_rel"], dtype=float)
    q_rad = q_rad / q_rad.mean()
    q_transport = np.asarray(data["q_transport_rel"], dtype=float)
    q_transport = q_transport / q_transport.mean()
    nwl = np.asarray(data["nwl_rel"], dtype=float)
    nwl = nwl / nwl.mean()
    d_min = float(np.min(data["d_m"]))

    # FWHM choices are resolution probes, not model claims. The wall standoff
    # is the key lower bound; the neutron width is included because it is the
    # surviving far-field channel's observed scale.
    widths = [
        ("raw withdrawn LOS", 0.0),
        ("FWHM = min standoff", d_min),
        ("FWHM = 2x min standoff", 2.0 * d_min),
        ("FWHM = neutron scale", 0.081),
    ]
    curves = [(name, fwhm, gaussian_periodic_smooth(q_rad, s, period, fwhm)) for name, fwhm in widths]

    summary = {
        "status": "diagnostic_only_not_a_promoted_v2_result",
        "source": str(args.input),
        "period_m": period,
        "sample_count": int(len(s)),
        "sample_spacing_m": ds,
        "minimum_standoff_m": d_min,
        "raw_radiative_fwhm_estimate_m": contiguous_fwhm_width(s, q_rad, period),
        "probe_curves": [
            {
                "label": name,
                "regularization_fwhm_m": fwhm,
                "radiative_peak_over_mean": float(curve.max()),
                "radiative_min_over_mean": float(curve.min()),
            }
            for name, fwhm, curve in curves
        ],
        "surviving_reference": {
            "cross_field_transport_peak_over_mean": float(q_transport.max()),
            "neutron_wall_load_min_over_mean": float(nwl.min()),
            "neutron_wall_load_max_over_mean": float(nwl.max()),
        },
        "acceptance_bar_before_v2_promotion": [
            "define near-field/far-field or analytic near-field treatment",
            "show radiative peak/mean is stable under wall-point and plasma-surface refinement",
            "emit a bounded OpenFOAM wall expression only after the convergence check passes",
        ],
    }

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_json = args.outdir / "radiative_nearfield_resolution_probe.json"
    out_png = args.outdir / "radiative_nearfield_resolution_probe.png"
    out_json.write_text(json.dumps(summary, indent=2))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12.2, 4.8), facecolor="white")
    x = 1000.0 * s
    for name, fwhm, curve in curves:
        lw = 2.6 if fwhm == 0 else 1.9
        ls = "--" if fwhm == 0 else "-"
        ax0.plot(x, curve, lw=lw, ls=ls, label=f"{name}: pk/mean {curve.max():.2f}")
    ax0.set_title("Withdrawn radiative LOS spike under length-scale probes")
    ax0.set_xlabel("duct arc length s [mm]")
    ax0.set_ylabel("radiative load / mean [-]")
    ax0.grid(alpha=0.28)
    ax0.legend(fontsize=8.2)

    labels = [name.replace("FWHM = ", "") for name, _, _ in curves]
    vals = [float(curve.max()) for _, _, curve in curves]
    colors = ["#b23b2b", "#5d8cc1", "#4f9d69", "#7b6fb0"]
    ax1.bar(np.arange(len(vals)), vals, color=colors)
    ax1.axhline(q_transport.max(), color="0.25", lw=1.2, ls=":", label="cross-field peak/mean 5.58")
    ax1.set_xticks(np.arange(len(vals)), labels, rotation=18, ha="right")
    ax1.set_ylabel("radiative peak/mean [-]")
    ax1.set_title("Regularization sensitivity, not a corrected result")
    ax1.grid(axis="y", alpha=0.28)
    ax1.legend(fontsize=8.4)

    fig.text(
        0.5,
        0.02,
        "Diagnostic only: v2 radiative/total peaking and hotspot extrapolation remain withdrawn until a converged near-field/far-field treatment exists.",
        ha="center",
        fontsize=9,
        color="0.25",
    )
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(out_png, dpi=150)

    # Live-tree runs mirror into the advisor packet. Packaged/flattened runs do
    # so only when the caller supplies --advisor-out.
    if args.advisor_out is not None:
        args.advisor_out.mkdir(parents=True, exist_ok=True)
        (args.advisor_out / "100_radiative_nearfield_resolution_probe.json").write_text(
            out_json.read_text()
        )
        import shutil

        shutil.copy2(
            out_png, args.advisor_out / "100_radiative_nearfield_resolution_probe.png"
        )
    print(json.dumps(summary, indent=2))
    print(f"wrote {out_png}")
    if args.advisor_out is not None:
        print(f"wrote {args.advisor_out / '100_radiative_nearfield_resolution_probe.png'}")


if __name__ == "__main__":
    main()
