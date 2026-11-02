#!/usr/bin/env python3
"""Stream-parse GPL17586 probe annotations and GSE110359 sample labels."""

from __future__ import annotations

import csv
import gzip
import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = Path("/Volumes/wang/paper15/ais_mesenchymal_axis/raw/gse110359")
SOFT = RAW / "GSE110359_family.soft.gz"
MATRIX = RAW / "GSE110359_series_matrix.txt.gz"
OUT = ROOT / "metadata/gse110359"


def unquote_fields(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))


def parse_samples() -> None:
    wanted: dict[str, list[str]] = {}
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith(("!Sample_title", "!Sample_geo_accession", "!Sample_source_name_ch1")):
                fields = unquote_fields(line)
                wanted[fields[0]] = fields[1:]
            elif line.startswith("!Sample_characteristics_ch1") and "disease state:" in line:
                fields = unquote_fields(line)
                wanted["disease_state"] = [value.split(":", 1)[1].strip() for value in fields[1:]]

    accessions = wanted["!Sample_geo_accession"]
    rows = []
    for index, accession in enumerate(accessions):
        rows.append(
            {
                "sample_id": accession,
                "title": wanted["!Sample_title"][index],
                "source": wanted["!Sample_source_name_ch1"][index],
                "disease_state": wanted["disease_state"][index].lower(),
            }
        )
    with (OUT / "sample_metadata.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def extract_gene_symbols(assignment: str) -> list[str]:
    symbols = []
    for record in assignment.split(" /// "):
        fields = record.split(" // ")
        if len(fields) >= 2:
            symbol = fields[1].strip()
            if symbol and symbol != "---" and symbol not in symbols:
                symbols.append(symbol)
    return symbols


def parse_platform() -> None:
    rows = []
    in_table = False
    header: list[str] | None = None
    with gzip.open(SOFT, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            record = dict(zip(header, fields))
            symbols = extract_gene_symbols(record.get("gene_assignment", ""))
            rows.append(
                {
                    "probe_id": record["ID"],
                    "gene_symbol": symbols[0] if len(symbols) == 1 else "",
                    "n_unique_gene_symbols": len(symbols),
                    "mapping_status": "unique" if len(symbols) == 1 else ("unmapped" if not symbols else "ambiguous"),
                    "category": record.get("category", ""),
                    "locus_type": record.get("locus type", ""),
                }
            )
    with (OUT / "probe_annotation.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    global RAW, SOFT, MATRIX, OUT
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=RAW)
    parser.add_argument("--out-dir", type=Path, default=OUT)
    args = parser.parse_args()
    RAW = args.raw_dir.expanduser().resolve()
    SOFT = RAW / "GSE110359_family.soft.gz"
    MATRIX = RAW / "GSE110359_series_matrix.txt.gz"
    OUT = args.out_dir.expanduser().resolve()
    for path in (SOFT, MATRIX):
        if not path.exists():
            raise FileNotFoundError(path)
    OUT.mkdir(parents=True, exist_ok=True)
    parse_samples()
    parse_platform()
    print("Parsed GSE110359 metadata and GPL17586 annotation")


if __name__ == "__main__":
    main()
