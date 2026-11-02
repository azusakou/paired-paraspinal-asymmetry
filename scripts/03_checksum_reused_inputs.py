#!/usr/bin/env python3
"""Record checksums for read-only public-data inputs reused from other projects."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bulk-counts", type=Path,
        default=Path("/Volumes/wang/paper14/ais_paraspinal_ecology/raw/bulk/GSE254300_readcount.xlsx"),
    )
    parser.add_argument(
        "--sc-object", type=Path,
        default=Path("/Volumes/wang/paper15/ais_mesenchymal_axis/derived/prjna722100/qc_annotation/prjna722100_qc_annotated_sce.rds"),
    )
    parser.add_argument(
        "--gse110359-dir", type=Path,
        default=Path("/Volumes/wang/paper15/ais_mesenchymal_axis/raw/gse110359"),
    )
    parser.add_argument("--output", type=Path, default=ROOT / "metadata/reused_input_manifest.tsv")
    args = parser.parse_args()
    inputs = {
        "GSE254300_readcount.xlsx": ("GSE254300", args.bulk_counts),
        "prjna722100_qc_annotated_sce.rds": ("PRJNA722100", args.sc_object),
        "GSE110359_series_matrix.txt.gz": ("GSE110359", args.gse110359_dir / "GSE110359_series_matrix.txt.gz"),
        "GSE110359_family.soft.gz": ("GSE110359", args.gse110359_dir / "GSE110359_family.soft.gz"),
    }
    rows = []
    for artifact, (dataset, path) in inputs.items():
        path = path.expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        rows.append(
            {
                "dataset": dataset,
                "artifact": artifact,
                "resolved_path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Recorded {len(rows)} reused inputs in {args.output}")


if __name__ == "__main__":
    main()
