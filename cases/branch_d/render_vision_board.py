#!/usr/bin/env python3
"""One-page VISION board to hook senior-design teammates. Leads with the two 3D results -- coolant FLOW
and first-wall HEAT (the step the ARC/Leffler deck listed as future work) -- then credibility (DNS) and
the pitch. Pure assembly of already-produced, honest artifacts."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg

plt.rcParams.update({"font.family": "DejaVu Sans"})
BG = "#0a0a0f"; FG = "#e8e8ee"; ACC = "#ff7a3d"; ACC2 = "#4dd0e1"

fig = plt.figure(figsize=(13, 17.5), facecolor=BG)
gs = fig.add_gridspec(4, 2, height_ratios=[0.52, 1.25, 0.95, 1.05], hspace=0.16, wspace=0.10,
                      left=0.045, right=0.955, top=0.97, bottom=0.03)

# --- title band ---
axt = fig.add_subplot(gs[0, :]); axt.axis("off")
axt.text(0.5, 0.80, "Simulating the heart of a fusion reactor", ha="center", va="center",
         color=FG, fontsize=30, fontweight="bold")
axt.text(0.5, 0.46, "Plasma -> neutrons -> coolant -> tritium: one pipeline, one honest slice",
         ha="center", va="center", color=ACC, fontsize=17)
axt.text(0.5, 0.12, "FLiBe-cooled stellarator first wall: the flow, heat, and breeding the ARC deck left as 'future work'",
         ha="center", va="center", color="#9a9aa8", fontsize=13, style="italic")

# --- hero 1: coolant flow cutaway ---
ax1 = fig.add_subplot(gs[1, :]); ax1.axis("off")
ax1.imshow(mpimg.imread("figs/vision_hero_coolant_slice.png"))

# --- hero 2: first-wall heat map (the differentiator) ---
ax2 = fig.add_subplot(gs[2, :]); ax2.axis("off")
ax2.imshow(mpimg.imread("figs/vision_firstwall_temp_plasma.png"))

# --- bottom-left: validation vs DNS ---
axv = fig.add_subplot(gs[3, 0]); axv.imshow(mpimg.imread("figs/opt1_c9_kader_apriori_Pr5.png")); axv.axis("off")
axv.set_title("Validated against published DNS", color=ACC2, fontsize=13, fontweight="bold", pad=4)

# --- bottom-right: pitch ---
axp = fig.add_subplot(gs[3, 1]); axp.axis("off")
bullets = [
    ("WHAT IT IS", FG, 13, "bold"),
    ("  3D wall-resolved LES of liquid-metal coolant +", "#c7c7d2", 11.5, "normal"),
    ("  heat through a curved fusion first-wall channel", "#c7c7d2", 11.5, "normal"),
    ("", FG, 6, "normal"),
    ("WHY IT'S REAL, NOT A TOY", FG, 13, "bold"),
    ("  - plasma-shaped wall load from a W7-X", "#c7c7d2", 11.5, "normal"),
    ("    equilibrium (DESC): hotspot 2.1x uniform", "#c7c7d2", 11.5, "normal"),
    ("  - OpenMC heating + tritium bred, transported,", "#c7c7d2", 11.5, "normal"),
    ("    mass balance closes to 0.02%", "#c7c7d2", 11.5, "normal"),
    ("  - benchmarked vs Kawamura + Moser DNS", "#c7c7d2", 11.5, "normal"),
    ("", FG, 6, "normal"),
    ("WHAT YOU'D BUILD / LEARN", FG, 13, "bold"),
    ("  OpenFOAM - OpenMC - HPC (Grace) -", "#c7c7d2", 11.5, "normal"),
    ("  Python - fusion thermal-hydraulics", "#c7c7d2", 11.5, "normal"),
    ("", FG, 6, "normal"),
    ("THE REACH", ACC, 13, "bold"),
    ("  scale this slice -> a full reactor first-wall", ACC, 11.5, "normal"),
    ("  segment on a national supercomputer", ACC, 11.5, "normal"),
]
y = 0.98
for txt, col, sz, w in bullets:
    axp.text(0.0, y, txt, color=col, fontsize=sz, fontweight=w, va="top", ha="left", transform=axp.transAxes)
    y -= 0.060 if txt else 0.030

fig.savefig("figs/vision_board.png", facecolor=BG, dpi=110, bbox_inches="tight")
print("wrote figs/vision_board.png")
