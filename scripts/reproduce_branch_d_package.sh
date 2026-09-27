#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib}"

echo "Rebuilding Branch D post-processing package from retained case data..."

python3 scripts/plot_branch_d_heat_flux_sensitivity.py
python3 scripts/plot_branch_d_fidelity_ladder_panel.py
python3 scripts/plot_branch_d_v27_mean_and_tke.py
python3 scripts/plot_branch_d_nearwall_flushing_law.py
python3 scripts/plot_branch_d_dimensionless_metrics.py
python3 scripts/plot_branch_d_v27_nearwall_recirculation.py
python3 scripts/plot_branch_d_bfs_reattachment_proxy.py

echo "Branch D package rebuilt."
echo "Key outputs:"
echo "  results/branch_d/"
echo "  results/report_figures/19_branch_d_heat_flux_sensitivity.png"
echo "  results/report_figures/21_branch_d_fidelity_ladder_panel.png"
echo "  results/report_figures/22_v27_rans_vs_les_mean_velocity_midplane.png"
echo "  results/report_figures/23_v27_resolved_tke_temperature_fluctuation_midplane.png"
echo "  results/report_figures/24_branch_d_nearwall_flushing_law.png"
echo "  results/report_figures/25_branch_d_dimensionless_metrics.png"
echo "  results/report_figures/26_v27_rans_vs_les_nearwall_recirculation_profile.png"
echo "  results/report_figures/27_v30_bfs_reattachment_proxy.png"
