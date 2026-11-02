#!/usr/bin/env python3
"""Build the manuscript Source Data bundle from analysis-native tables."""

from __future__ import annotations

import csv
import hashlib
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source_data"

FILES = {
    "Table1_dataset_registry.tsv": "metadata/dataset_registry.tsv",
    "Figure1_gene_effects.tsv": "results/integration/orthogonal_axes/gene_evidence_axis_map.tsv",
    "Figure1_model_comparison.tsv": "results/manuscript_figures_nature/Figure1_model_comparison.tsv",
    "Model_sensitivity_summary.tsv": "results/tissue/model_dependence/model_summary.tsv",
    "Equivalence_model_comparison.tsv": "results/tissue/model_dependence/equivalence_model_comparison.tsv",
    "Figure2a_gene_equivalence.tsv.gz": "results/tissue/gse254300_bulk/cross_etiology_equivalence/gene_cross_etiology_equivalence.tsv.gz",
    "Figure2b_equivalence_margins.tsv": "results/tissue/gse254300_bulk/cross_etiology_equivalence/equivalence_margin_sensitivity.tsv",
    "Figure2c_leave_one_donor_out.tsv": "results/tissue/gse254300_bulk/leave_one_donor_out_correlations.tsv",
    "Figure3_camera_pathways.tsv": "results/integration/orthogonal_axes/source_data_figure3.tsv",
    "Figure4a_composition_sensitivity.tsv": "results/integration/composition_precision_sensitivity/source_data_figure4a.tsv",
    "Figure4b_celltype_pathway_scores.tsv": "results/integration/orthogonal_axes/celltype_pathways/donor_paired_celltype_pathway_effects.tsv",
    "Figure5_origin_class_counts.tsv": "results/integration/origin_annotation/origin_evidence_class_counts.tsv",
    "Figure5_source_threshold_sensitivity.tsv": "results/integration/origin_annotation/source_priority_threshold_sensitivity.tsv",
    "Figure5a_evidence_funnel.tsv": "results/manuscript_figures_nature/Figure5_funnel.tsv",
    "Figure5c_bridge_candidates.tsv": "results/integration/orthogonal_axes/shared_adaptation_plasma_bridge_candidates.tsv",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for destination, source in FILES.items():
        source_path = ROOT / source
        if not source_path.exists():
            raise FileNotFoundError(source_path)
        destination_path = OUT / destination
        shutil.copy2(source_path, destination_path)
        rows.append(
            {
                "source_data_file": destination,
                "analysis_source": source,
                "bytes": destination_path.stat().st_size,
                "sha256": sha256(destination_path),
            }
        )
    with (OUT / "manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Built {len(rows)} Source Data files in {OUT}")


if __name__ == "__main__":
    main()
