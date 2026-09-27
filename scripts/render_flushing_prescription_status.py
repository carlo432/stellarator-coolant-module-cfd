#!/usr/bin/env python3
"""Replace the old local-flushing prescription figure with its tested status."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "advisor_packet" / "96_flushing_remedy.png"


def main() -> None:
    predicted = 1.221
    measured = 1.033

    fig, (ax, note) = plt.subplots(
        1,
        2,
        figsize=(12.4, 5.1),
        facecolor="white",
        gridspec_kw={"width_ratios": [1.05, 1.35]},
    )
    bars = ax.bar(
        [0, 1],
        [predicted, measured],
        color=["#9a6b2f", "#1d7059"],
        width=0.58,
    )
    ax.axhline(1.0, color="0.45", ls="--", lw=1.2)
    ax.set_xticks([0, 1], ["developed-law\nprediction", "34 mm throat\nCFD result"])
    ax.set_ylabel("local superheat reduction factor")
    ax.set_ylim(0.98, 1.26)
    ax.set_title("Pointwise prescription under-delivers in a real throat", loc="left")
    ax.grid(axis="y", alpha=0.25)
    for bar, value in zip(bars, [predicted, measured]):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.006, f"{value:.3f}x", ha="center", weight="bold")

    note.axis("off")
    note.set_title("Prescription status: local implementation retired", loc="left", fontsize=13)
    note.text(
        0,
        0.88,
        "The developed CHT sweep established dT ∝ U⁻⁰·⁷⁰.\n"
        "Applied pointwise, it predicted 18% local cooling\n"
        "from a 1.33x speed increase.\n\n"
        "The direct throat case delivered about 1.30x near-wall\n"
        "acceleration but only 3.3% local cooling. P99 film\n"
        "temperature remained 1898 K, above 1703 K boiling.\n\n"
        "RETAIN: developed-flow law.\n"
        "RETIRE: 34 mm local throat as a boiling remedy.\n"
        "NEXT: whole-duct/long-section speed or channel resize.",
        va="top",
        fontsize=11.2,
        linespacing=1.35,
    )
    fig.suptitle("Flushing prescription after direct CFD test", fontsize=16, weight="bold")
    fig.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
