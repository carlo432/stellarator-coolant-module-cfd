#!/usr/bin/env python3
"""Check which Kawamura DNS digitization CSVs have real data rows."""
from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "references" / "digitization" / "kawamura_1998"
MANIFEST = WORK / "target_manifest.csv"


def count_data_rows(path: Path) -> int:
    if not path.exists():
        return -1
    with path.open(newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        return 0
    return max(len(rows) - 1, 0)


def main() -> None:
    with MANIFEST.open(newline="") as f:
        targets = list(csv.DictReader(f))

    print("figure_id,priority,digitized_csv,data_rows,status")
    missing_required = 0
    for row in targets:
        path = WORK / row["digitized_csv"]
        rows = count_data_rows(path)
        if rows < 0:
            status = "missing_file"
        elif rows == 0:
            status = "empty"
        else:
            status = "has_data"
        if row["priority"] == "required" and status != "has_data":
            missing_required += 1
        print(f"{row['figure_id']},{row['priority']},{row['digitized_csv']},{rows},{status}")

    print(f"required_targets_missing_data={missing_required}")
    if missing_required:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
