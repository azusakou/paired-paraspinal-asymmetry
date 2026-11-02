#!/usr/bin/env python3
"""Fail-fast validation of the conserved-paraspinal-adaptation pivot."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def rows(path: str):
    target = ROOT / path
    opener = gzip.open if target.suffix == ".gz" else open
    with opener(target, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def metric_map(path: str):
    return {row["metric"]: float(row["value"]) for row in rows(path)}


def require(condition: bool, message: str):
    if not condition:
        raise AssertionError(message)
    print(f"PASS\t{message}")


def require_current_output(output: str, producer: str):
    output_path = ROOT / output
    producer_path = ROOT / producer
    require(output_path.exists(), f"output exists: {output}")
    require(
        output_path.stat().st_mtime >= producer_path.stat().st_mtime,
        f"output is not older than producer: {output}",
    )


axis = metric_map("results/integration/orthogonal_axes/analysis_summary.tsv")
require(axis["shared_adaptation_genes_fdr05"] == 743, "743 shared-effect genes pass FDR 0.05")
require(axis["shared_adaptation_same_direction"] == 742, "742 shared hits have concordant AIS/CS direction")
require(axis["interaction_genes_fdr10"] == 0, "no disease-by-side interaction passes FDR 0.10")
require(axis["msc_genes_fdr05"] == 0 and axis["msc_genes_fdr10"] == 4, "BM-MSC gene-level signal remains weak")
require(axis["shared_pathways_camera_fdr05"] == 25, "25 shared-side CAMERA pathways pass FDR 0.05")
require(axis["msc_pathways_camera_fdr05"] == 70, "70 BM-MSC CAMERA pathways pass FDR 0.05")
require(axis["joint_concordant_pathways_camera_fdr05"] == 0, "no jointly significant concordant CAMERA pathway")
require(axis["joint_opposed_pathways_camera_fdr05"] == 1, "one jointly significant opposed CAMERA pathway")

etiology_pathways = metric_map("results/tissue/gse254300_bulk/camera_cross_etiology_pathway_summary.tsv")
require(
    etiology_pathways["AIS_camera_fdr05"] == 29
    and etiology_pathways["CS_camera_fdr05"] == 31,
    "separate CAMERA analyses identify 29 AIS and 31 congenital-scoliosis pathways",
)
require(
    etiology_pathways["joint_AIS_CS_camera_fdr05"] == 9
    and etiology_pathways["joint_AIS_CS_camera_fdr05_same_direction"] == 8,
    "nine pathways are FDR-significant in both etiologies and eight share direction",
)
require(
    etiology_pathways["mean_side_hits_nominal_p05_in_both"] == 25
    and etiology_pathways["mean_side_hits_same_direction_in_both"] == 25,
    "all 25 mean-side pathways are nominally supported with the same direction in both strata",
)

loo = metric_map("results/tissue/gse254300_bulk/leave_one_donor_out_summary.tsv")
require(loo["full_shared_fdr05_all_loo_same_direction"] == 743, "all 743 shared hits retain direction in every leave-one-donor-out fit")
require(loo["minimum_loo_spearman"] > 0.95, "minimum leave-one-donor-out effect correlation exceeds 0.95")

perm = metric_map("results/tissue/gse254300_bulk/permutation_sensitivity/permutation_summary.tsv")
require(perm["minimum_two_sided_sign_flip_p"] == 2 / 1024, "exact sign-flip minimum P value matches 2/1024 resolution")
require(perm["edgeR_shared_fdr05_exact_p_le_0.01"] == 678, "678 edgeR hits have exact P <= 0.01")
require(perm["edgeR_shared_fdr05_exact_p_le_0.05"] == 736, "736 edgeR hits have exact P <= 0.05")
require(perm["effect_spearman"] > 0.97, "edgeR and paired-logCPM effects correlate above 0.97")
require(perm["edgeR_interaction_fdr10"] == 0 and perm["exact_group_permutation_interaction_fdr10"] == 0, "interaction null result is method-robust")

bridge = rows("results/integration/orthogonal_axes/shared_adaptation_plasma_bridge_candidates.tsv")
bridge_genes = {row["gene_symbol"] for row in bridge}
require(bridge_genes == {"ECM1", "LYVE1", "PI16"}, "bridge set is exactly ECM1, LYVE1 and PI16")

exact = {row["gene_symbol"]: row for row in rows("results/tissue/gse254300_bulk/permutation_sensitivity/exact_permutation_results.tsv")}
loo_gene = {row["gene_symbol"]: row for row in rows("results/tissue/gse254300_bulk/leave_one_donor_out_shared_effects.tsv.gz")}
for gene in sorted(bridge_genes):
    require(float(exact[gene]["sign_flip_PValue"]) <= 4 / 1024, f"{gene} exact paired P <= 4/1024")
    require(loo_gene[gene]["all_loo_same_direction"] == "TRUE", f"{gene} retains direction in all leave-one-donor-out fits")

equivalence = metric_map("results/tissue/gse254300_bulk/cross_etiology_equivalence/equivalence_summary.tsv")
require(equivalence["strictly_conserved_margin_0p50"] == 163, "163 shared genes meet FDR-controlled +/-0.50-log2 equivalence")
require(0.21 < equivalence["strictly_conserved_fraction_of_shared"] < 0.23, "strict equivalence comprises about 22% of shared hits")
margin_rows = rows("results/tissue/gse254300_bulk/cross_etiology_equivalence/equivalence_margin_sensitivity.tsv")
margin_counts = {float(row["equivalence_margin_log2"]): int(row["shared_hits_equivalent_fdr05"]) for row in margin_rows}
require(margin_counts == {0.25: 0, 0.5: 163, 1.0: 590}, "equivalence sensitivity is 0/163/590 at +/-0.25/0.50/1.00 log2 margins")

composition = metric_map("results/integration/composition_precision_sensitivity/analysis_summary.tsv")
require(composition["original_camera_pathways_fdr05"] == 25, "composition sensitivity reproduces 25 primary CAMERA pathways")
require(composition["retained_after_excluding_specificity_ge_0p70_fdr05"] == 21, "21 pathways survive exclusion of specificity >=0.70 genes")
require(composition["retained_after_excluding_specificity_ge_0p50_fdr05"] == 21, "21 pathways survive exclusion of specificity >=0.50 genes")
require(composition["precision_adjusted_pathways_fdr05"] == 0, "strict-equivalence subset has no precision-matched pathway enrichment")
require_current_output(
    "results/integration/orthogonal_axes/orthogonal_axis_pathways_camera.pdf",
    "scripts/22_analyze_orthogonal_evidence_axes.R",
)
require_current_output(
    "results/integration/composition_precision_sensitivity/celltype_specificity_pathway_sensitivity.pdf",
    "scripts/27_test_composition_and_precision_sensitivity.R",
)
require_current_output(
    "results/integration/origin_annotation/origin_evidence_class_counts.pdf",
    "scripts/15_integrate_external_origin_annotations.R",
)

power_targets = rows("results/design/prospective_validation/etiology_equivalence_sample_targets.tsv")
target_map = {
    (float(row["standardized_equivalence_margin"]), float(row["target_power"])): int(row["complete_donors_per_etiology"])
    for row in power_targets
}
require(target_map[(0.5, 0.8)] == 70 and target_map[(0.5, 0.9)] == 88, "strict +/-0.50-SD equivalence requires 70/88 complete donors per etiology for 80/90% power")
require(target_map[(0.75, 0.8)] == 32 and target_map[(0.75, 0.9)] == 40, "wider +/-0.75-SD equivalence requires 32/40 complete donors per etiology for 80/90% power")

gates = rows("metadata/prediction_gate.tsv")
require(all(row["status"] == "HOLD" for row in gates), "all prediction gates remain HOLD")

print("VALIDATION_COMPLETE\tconserved-paraspinal-adaptation pivot")
