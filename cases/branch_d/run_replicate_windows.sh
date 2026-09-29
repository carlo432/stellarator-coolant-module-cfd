#!/usr/bin/env bash
# Replicate averaging windows on the plasma baseline = the reference every MHD / k2 /
# sensitivity case is differenced against. Two consecutive independent 0.7 s windows on a
# fully developed field give a MEASURED window-to-window scatter for sup99 and sigma_hot,
# replacing the assumed "+-2-3%".
set -e
cd "$(dirname "$0")"
until grep -q DONE geom_nowave_drift_ctl2/run.out 2>/dev/null; do sleep 15; done
python3 make_drift_control.py geom_pipeline_plasma      geom_plasma_rep1 4.3
bash geom_plasma_rep1/run_pipeline.sh   > geom_plasma_rep1/run.out 2>&1
python3 make_drift_control.py geom_plasma_rep1          geom_plasma_rep2 5.0
bash geom_plasma_rep2/run_pipeline.sh   > geom_plasma_rep2/run.out 2>&1
echo ALLDONE
