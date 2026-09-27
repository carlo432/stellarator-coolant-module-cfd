#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CASE="$ROOT/cases/branch_d/c9_thermal_cert_run"
LOG="$CASE/log.pimpleFoam"
SENTINEL="$CASE/C9.done"
END_TIME="1.2"

print_status() {
  local now status latest progress mtime size last_write
  now="$(date '+%Y-%m-%d %H:%M:%S %Z')"

  if [[ -f "$SENTINEL" ]]; then
    status="DONE sentinel present"
  elif [[ -f "$LOG" ]] && tail -n 120 "$LOG" | rg -q '^End$|Finalising parallel run'; then
    status="SOLVER REACHED END; reconstruction or C9.done sentinel still pending"
  else
    status="RUNNING or finishing; no C9.done sentinel yet"
  fi

  if [[ -f "$LOG" ]]; then
    latest="$(rg '^Time = ' "$LOG" | tail -n 1 | awk '{print $3}')"
    mtime="$(stat -c '%y' "$LOG")"
    size="$(stat -c '%s bytes' "$LOG")"
    progress="$(awk -v t="${latest:-0}" -v e="$END_TIME" 'BEGIN { if (e > 0) printf "%.1f%%", 100*t/e; else print "n/a" }')"
  else
    latest="n/a"
    mtime="missing log"
    size="missing log"
    progress="n/a"
  fi

  last_write="$(find "$CASE" -type f -printf '%T@ %TY-%Tm-%Td %TH:%TM:%TS %p\n' 2>/dev/null | sort -n | tail -n 1 | cut -d' ' -f2-)"

  cat <<EOF
# C9 Live Status

Checked: $now

- status: $status
- latest logged solver time: $latest / $END_TIME
- approximate progress to configured endTime: $progress
- log mtime: $mtime
- log size: $size
- newest file write under case: ${last_write:-n/a}

Case:

\`$CASE\`

To watch continuously:

\`\`\`bash
scripts/watch_c9_status.sh --loop
\`\`\`
EOF
}

if [[ "${1:-}" == "--loop" ]]; then
  while true; do
    clear
    print_status
    sleep "${WATCH_SECONDS:-10}"
  done
else
  print_status
fi
