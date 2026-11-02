#!/usr/bin/env python3
"""Fail closed when Phase-0 evidence is incomplete or internally inconsistent."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.project_root.resolve()

    required = [
        root / "metadata" / "dataset_registry.tsv",
        root / "metadata" / "data_dictionary.tsv",
        root / "metadata" / "prediction_gate.tsv",
        root / "metadata" / "published_overlap.tsv",
        root / "metadata" / "pxd024246_published_methods.tsv",
        root / "metadata" / "retrieval_log.tsv",
        root / "reports" / "phase0_data_availability_audit.md",
        root / "reports" / "phase0_fair_audit.md",
        root / "results" / "audit" / "pxd024246" / "pxd024246_parse_audit.json",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("Missing Phase-0 artifacts:\n" + "\n".join(missing))

    gates = read_tsv(root / "metadata" / "prediction_gate.tsv")
    open_gates = [row for row in gates if row["status"] != "PASS"]

    with (root / "results" / "audit" / "pxd024246" / "pxd024246_parse_audit.json").open(
        encoding="utf-8"
    ) as handle:
        audit = json.load(handle)

    expected = {
        "abundance_columns": 84,
        "px_raw_file_count": 84,
        "imputation_applied": False,
    }
    mismatches = {
        key: {"expected": value, "observed": audit.get(key)}
        for key, value in expected.items()
        if audit.get(key) != value
    }
    if mismatches:
        raise SystemExit("Phase-0 audit mismatch: " + json.dumps(mismatches, sort_keys=True))

    report = {
        "phase0_artifacts_present": True,
        "prediction_gate": "HOLD" if open_gates else "PASS",
        "open_gate_count": len(open_gates),
        "open_gates": [row["gate"] for row in open_gates],
        "abundance_columns": audit["abundance_columns"],
        "raw_files_listed": audit["px_raw_file_count"],
        "unmatched_abundance_columns": audit["abundance_columns_missing_from_px_raw_list"],
        "unmatched_raw_files": audit["px_raw_files_missing_from_abundance_columns"],
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
