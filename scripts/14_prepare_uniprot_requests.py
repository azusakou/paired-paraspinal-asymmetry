#!/usr/bin/env python3
"""Prepare reproducible UniProt REST request payloads for PXD024246 accessions."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "results/integration/tissue_support/plasma_protein_tissue_support_matrix.tsv"
OUTDIR = ROOT / "metadata/uniprot_requests"
EXTERNAL_ROOT = Path(
    os.environ.get(
        "PAPER16_EXTERNAL_ROOT", "/Volumes/wang/paper16/ais_sc_plasma/raw/external"
    )
)
RAW_DIR = EXTERNAL_ROOT / "UniProt_2026-07-22"
CHUNK_SIZE = 40


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    with MATRIX.open(newline="", encoding="utf-8") as handle:
        accessions = sorted(
            {
                row["accession"].strip()
                for row in csv.DictReader(handle, delimiter="\t")
                if row.get("accession", "").strip()
            }
        )

    manifest = []
    for index, start in enumerate(range(0, len(accessions), CHUNK_SIZE), start=1):
        chunk = accessions[start : start + CHUNK_SIZE]
        payload_path = OUTDIR / f"uniprot_chunk_{index:02d}.json"
        raw_path = RAW_DIR / f"uniprot_chunk_{index:02d}.tsv"
        payload = {
            "base_url": "https://rest.uniprot.org",
            "path": "uniprotkb/search",
            "params": {
                "query": " OR ".join(f"accession:{accession}" for accession in chunk),
                "format": "tsv",
                "size": 500,
                "fields": (
                    "accession,id,gene_names,protein_name,organism_name,reviewed,"
                    "cc_subcellular_location,ft_signal,cc_tissue_specificity,"
                    "cc_function,go_c"
                ),
            },
            "response_format": "text",
            "timeout_sec": 60,
            "save_raw": True,
            "raw_output_path": str(raw_path),
            "max_items": 5,
            "max_depth": 3,
        }
        payload_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        manifest.append(
            {
                "chunk": index,
                "n_accessions": len(chunk),
                "payload": str(payload_path.relative_to(ROOT)),
                "raw_output": str(raw_path),
            }
        )

    with (OUTDIR / "manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Prepared {len(manifest)} requests for {len(accessions)} unique accessions")


if __name__ == "__main__":
    main()
