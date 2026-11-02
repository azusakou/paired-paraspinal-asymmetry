#!/usr/bin/env python3
"""Parse the PXD024246 Proteome Discoverer SpreadsheetML export.

The script intentionally performs no transformation, normalization, filtering,
or imputation. Blank SpreadsheetML cells remain missing values.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import xml.etree.ElementTree as ET
from pathlib import Path

SS = "urn:schemas-microsoft-com:office:spreadsheet"
NS = {"ss": SS}
INDEX = f"{{{SS}}}Index"


def parse_row(row: ET.Element, width: int) -> list[str | None]:
    values: list[str | None] = [None] * width
    position = 0
    for cell in row.findall("ss:Cell", NS):
        explicit = cell.attrib.get(INDEX)
        if explicit is not None:
            position = int(explicit) - 1
        data = cell.find("ss:Data", NS)
        values[position] = None if data is None else data.text
        position += 1
    return values


def numeric(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def raw_names_from_px(px_xml: Path) -> set[str]:
    root = ET.parse(px_xml).getroot()
    names = set()
    for element in root.iter("DatasetFile"):
        name = element.attrib.get("name", "")
        if name.lower().endswith(".raw"):
            names.add(Path(name).stem)
    return names


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--proteins-xml", required=True, type=Path)
    parser.add_argument("--px-xml", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    tree = ET.parse(args.proteins_xml)
    root = tree.getroot()
    table = root.find(".//ss:Worksheet[@ss:Name='Proteins']/ss:Table", NS)
    if table is None:
        raise RuntimeError("Proteins worksheet not found")
    width = int(table.attrib[f"{{{SS}}}ExpandedColumnCount"])
    rows = table.findall("ss:Row", NS)
    parsed = [parse_row(row, width) for row in rows]
    header = ["" if value is None else value for value in parsed[0]]
    fixed_width = 18
    fixed = header[:fixed_width]
    sample_ids = header[fixed_width:]
    if not all(sample_ids):
        raise RuntimeError("One or more abundance columns have no sample ID")

    feature_rows: list[list[object]] = []
    matrix_rows: list[list[object]] = []
    missing_by_sample = [0] * len(sample_ids)
    observed_by_sample = [0] * len(sample_ids)
    master_count = 0
    high_confidence_count = 0

    for row in parsed[1:]:
        meta = row[:fixed_width]
        abundance_text = row[fixed_width:]
        abundances = [numeric(value) for value in abundance_text]
        if meta[2] == "Master Protein":
            master_count += 1
        if meta[1] == "High":
            high_confidence_count += 1
        for i, value in enumerate(abundances):
            if value is None:
                missing_by_sample[i] += 1
            else:
                observed_by_sample[i] += 1
        feature_rows.append([
            meta[3], meta[4], meta[1], meta[2], meta[5], meta[7], meta[8],
            meta[9], meta[10], meta[11], meta[13], meta[14],
            sum(value is None for value in abundances),
            sum(value is not None for value in abundances),
        ])
        matrix_rows.append([meta[3]] + ["" if v is None else format(v, ".17g") for v in abundances])

    write_tsv(
        args.outdir / "protein_features.tsv",
        ["accession", "description", "fdr_confidence", "master_status", "q_value",
         "coverage_percent", "n_peptides", "n_psms", "n_unique_peptides",
         "n_protein_groups", "mw_kda", "calc_pi", "n_missing", "n_observed"],
        feature_rows,
    )
    write_tsv(args.outdir / "abundance_matrix_raw_export.tsv", ["accession"] + sample_ids, matrix_rows)

    px_raw = raw_names_from_px(args.px_xml)
    sample_set = set(sample_ids)
    sample_rows = []
    n_proteins = len(matrix_rows)
    for sample, missing, observed in zip(sample_ids, missing_by_sample, observed_by_sample):
        sample_rows.append([
            sample, f"{sample}.raw" if sample in px_raw else "",
            observed, missing, missing / n_proteins,
        ])
    write_tsv(
        args.outdir / "sample_missingness.tsv",
        ["abundance_column", "matched_raw_file", "n_observed", "n_missing", "missing_fraction"],
        sample_rows,
    )

    fractions = [missing / n_proteins for missing in missing_by_sample]
    audit = {
        "source_file": str(args.proteins_xml),
        "spreadsheet_columns": width,
        "annotation_columns": fixed_width,
        "abundance_columns": len(sample_ids),
        "protein_rows_all": n_proteins,
        "master_protein_rows": master_count,
        "high_confidence_rows": high_confidence_count,
        "px_raw_file_count": len(px_raw),
        "abundance_columns_missing_from_px_raw_list": sorted(sample_set - px_raw, key=lambda x: int(x)),
        "px_raw_files_missing_from_abundance_columns": sorted(px_raw - sample_set, key=lambda x: int(x)),
        "sample_missing_fraction": {
            "min": min(fractions),
            "median": statistics.median(fractions),
            "max": max(fractions),
        },
        "transformations_applied": [],
        "imputation_applied": False,
    }
    (args.outdir / "pxd024246_parse_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

