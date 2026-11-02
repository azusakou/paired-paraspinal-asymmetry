#!/usr/bin/env python3
"""Validate Phase-1 tissue/source outputs and modelling guardrails."""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXTERNAL = Path(
    os.environ.get(
        "PAPER16_EXTERNAL_ROOT", "/Volumes/wang/paper16/ais_sc_plasma/raw/external"
    )
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def truth(value: str) -> bool:
    return value.upper() == "TRUE"


def main() -> None:
    errors: list[str] = []
    origin_path = ROOT / "results/integration/origin_annotation/plasma_protein_origin_evidence_matrix.tsv"
    counts_path = ROOT / "results/integration/origin_annotation/origin_evidence_class_counts.tsv"
    candidates_path = ROOT / "results/integration/origin_annotation/high_confidence_paraspinal_secretory_compatible.tsv"
    sensitivity_path = ROOT / "results/integration/origin_annotation/source_priority_threshold_sensitivity.tsv"
    for path in (origin_path, counts_path, candidates_path, sensitivity_path):
        if not path.exists():
            errors.append(f"Missing output: {path}")
    if errors:
        raise SystemExit("\n".join(errors))

    origin = read_tsv(origin_path)
    counts = read_tsv(counts_path)
    candidates = read_tsv(candidates_path)
    sensitivity = read_tsv(sensitivity_path)
    if len(origin) != 394:
        errors.append(f"Expected 394 protein rows, found {len(origin)}")
    if len({row["accession"] for row in origin}) != 394:
        errors.append("Origin matrix accessions are not unique")
    if sum(int(row["n_proteins"]) for row in counts) != len(origin):
        errors.append("Origin class counts do not sum to origin matrix rows")
    allowed_support = {
        "strong_celltype_localization", "moderate_celltype_localization",
        "broad_expression", "low_or_sparse_expression",
        "not_detected_in_sc_reference", "unmapped_protein_export", "NA", "",
    }
    unexpected_support = {
        row["paraspinal_expression_support"] for row in origin
        if row["paraspinal_expression_support"] not in allowed_support
    }
    if unexpected_support:
        errors.append(f"Unrecognized paraspinal support labels: {sorted(unexpected_support)}")
    for row in origin:
        if (
            row["paraspinal_expression_support"] in
            {"not_detected_in_sc_reference"}
            and not truth(row["competing_systemic_source"])
            and row["cell_origin_evidence_class"] != "3_plasma_detected_paraspinal_low_or_absent"
        ):
            errors.append(
                f"Low/absent paraspinal entry {row['accession']} escaped class 3"
            )
    scenarios = {row["scenario"] for row in sensitivity}
    if scenarios != {"permissive", "primary", "stringent"}:
        errors.append(f"Source-priority sensitivity scenarios incomplete: {sorted(scenarios)}")

    for row in candidates:
        if row["cell_origin_evidence_class"] != "1_paraspinal_secretory_compatible_not_proven":
            errors.append(f"Candidate {row['accession']} has wrong evidence class")
        if row["fdr_confidence"] != "High":
            errors.append(f"Candidate {row['accession']} is not a high-confidence PXD master protein")
        if truth(row["competing_systemic_source"]):
            errors.append(f"Candidate {row['accession']} retains prioritized systemic-source evidence")
        if not truth(row["classical_or_explicit_secretory_support"]):
            errors.append(f"Candidate {row['accession']} lacks secretory support")

    prediction_gate = read_tsv(ROOT / "metadata/prediction_gate.tsv")
    if any(row["status"] != "HOLD" for row in prediction_gate):
        errors.append("Prediction gate was opened without a documented Phase-1 label resolution")

    manifest = read_tsv(ROOT / "metadata/download_manifest.tsv")
    manifest_by_artifact = {row["artifact"]: row for row in manifest}
    hpa_files = [
        EXTERNAL / "HPA_v25.1/rna_tissue_consensus.tsv.zip",
        EXTERNAL / "HPA_v25.1/rna_immune_cell.tsv.zip",
    ]
    for path in hpa_files:
        if not path.exists():
            errors.append(f"Missing external input: {path}")
            continue
        record = manifest_by_artifact.get(path.name)
        if record is None:
            errors.append(f"No manifest record for {path.name}")
        elif sha256(path) != record["sha256"]:
            errors.append(f"Checksum mismatch for {path.name}")

    uniprot_files = sorted((EXTERNAL / "UniProt_2026-07-22").glob("uniprot_chunk_*.tsv"))
    if len(uniprot_files) != 10:
        errors.append(f"Expected 10 UniProt chunks, found {len(uniprot_files)}")
    returned_accessions: list[str] = []
    for path in uniprot_files:
        rows = read_tsv(path)
        returned_accessions.extend(row["Entry"] for row in rows)
        record = manifest_by_artifact.get(path.name)
        if record is None or sha256(path) != record["sha256"]:
            errors.append(f"Missing or mismatched manifest checksum for {path.name}")
    if len(returned_accessions) != 394 or len(set(returned_accessions)) != 394:
        errors.append("UniProt response is not a one-to-one set of 394 accessions")

    if errors:
        raise SystemExit("PHASE 1 VALIDATION FAILED\n- " + "\n- ".join(errors))
    print(
        "PHASE 1 VALIDATION PASSED: 394 proteins; "
        f"{len(candidates)} high-confidence paraspinal-secretory-compatible candidates; "
        "prediction remains HOLD"
    )


if __name__ == "__main__":
    main()
