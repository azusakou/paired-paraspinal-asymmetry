#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 160)

file_arg <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])
root <- normalizePath(file.path(dirname(file_arg), ".."), mustWork = TRUE)
matrix_path <- file.path(root, "results", "integration", "tissue_support", "plasma_protein_tissue_support_matrix.tsv")
external_root <- Sys.getenv("PAPER16_EXTERNAL_ROOT", "/Volumes/wang/paper16/ais_sc_plasma/raw/external")
hpa_dir <- file.path(external_root, "HPA_v25.1")
uniprot_dir <- file.path(external_root, "UniProt_2026-07-22")
out_dir <- file.path(root, "results", "integration", "origin_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  matrix_path,
  file.path(hpa_dir, "rna_tissue_consensus.tsv.zip"),
  file.path(hpa_dir, "rna_immune_cell.tsv.zip")
)
missing_inputs <- required[!file.exists(required)]
if (length(missing_inputs)) stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "))

plasma <- read.delim(matrix_path, check.names = FALSE)
target_ensembl <- unique(na.omit(plasma$ensembl_id))

read_zip_tsv <- function(zip_path, member) {
  extract_dir <- tempfile("paper16_hpa_")
  dir.create(extract_dir)
  on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)
  utils::unzip(zip_path, files = member, exdir = extract_dir)
  read.delim(file.path(extract_dir, member), check.names = FALSE)
}

tissue <- read_zip_tsv(file.path(hpa_dir, "rna_tissue_consensus.tsv.zip"), "rna_tissue_consensus.tsv")
tissue <- tissue[tissue$Gene %in% target_ensembl, , drop = FALSE]
names(tissue)[names(tissue) == "nTPM"] <- "ntpm"

immune <- read_zip_tsv(file.path(hpa_dir, "rna_immune_cell.tsv.zip"), "rna_immune_cell.tsv")
immune <- immune[immune$Gene %in% target_ensembl, , drop = FALSE]
names(immune)[names(immune) == "nTPM"] <- "ntpm"

summarise_tissue <- function(df) {
  split_rows <- split(df, df$Gene)
  do.call(rbind, lapply(split_rows, function(x) {
    values <- setNames(x$ntpm, x$Tissue)
    decreasing_order <- order(x$ntpm, decreasing = TRUE, na.last = TRUE)
    top <- x[decreasing_order[1], , drop = FALSE]
    second <- if (nrow(x) >= 2) x$ntpm[decreasing_order[2]] else NA_real_
    selected <- function(name) if (name %in% names(values)) unname(values[name]) else NA_real_
    data.frame(
      ensembl_id = x$Gene[1],
      hpa_gene_name = x$`Gene name`[1],
      hpa_tissue_max = top$Tissue[1],
      hpa_tissue_max_ntpm = top$ntpm[1],
      hpa_tissue_second_ntpm = second,
      hpa_liver_ntpm = selected("liver"),
      hpa_skeletal_muscle_ntpm = selected("skeletal muscle"),
      hpa_smooth_muscle_ntpm = selected("smooth muscle"),
      hpa_adipose_ntpm = selected("adipose tissue"),
      hpa_bone_marrow_ntpm = selected("bone marrow"),
      hpa_spleen_ntpm = selected("spleen"),
      hpa_lymph_node_ntpm = selected("lymph node"),
      stringsAsFactors = FALSE
    )
  }))
}

summarise_immune <- function(df) {
  split_rows <- split(df, df$Gene)
  do.call(rbind, lapply(split_rows, function(x) {
    x_non_pbmc <- x[x$`Immune cell` != "total PBMC", , drop = FALSE]
    top <- x_non_pbmc[which.max(x_non_pbmc$ntpm), , drop = FALSE]
    data.frame(
      ensembl_id = x$Gene[1],
      hpa_immune_max_cell = top$`Immune cell`[1],
      hpa_immune_max_ntpm = top$ntpm[1],
      hpa_total_pbmc_ntpm = if (any(x$`Immune cell` == "total PBMC")) x$ntpm[x$`Immune cell` == "total PBMC"][1] else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

tissue_summary <- summarise_tissue(tissue)
immune_summary <- summarise_immune(immune)

uniprot_files <- sort(list.files(uniprot_dir, pattern = "^uniprot_chunk_[0-9]+\\.tsv$", full.names = TRUE))
if (!length(uniprot_files)) stop("No UniProt chunk TSV files found in ", uniprot_dir)
uniprot <- do.call(rbind, lapply(uniprot_files, function(path) read.delim(path, check.names = FALSE, quote = "", comment.char = "")))
if (anyDuplicated(uniprot$Entry)) stop("Duplicated UniProt entries across chunks")

uniprot_small <- data.frame(
  accession = uniprot$Entry,
  uniprot_reviewed = uniprot$Reviewed,
  uniprot_primary_gene = sub(" .*", "", uniprot$`Gene Names`),
  uniprot_subcellular_location = uniprot$`Subcellular location [CC]`,
  uniprot_signal_peptide = uniprot$`Signal peptide`,
  uniprot_tissue_specificity = uniprot$`Tissue specificity`,
  uniprot_go_cellular_component = uniprot$`Gene Ontology (cellular component)`,
  stringsAsFactors = FALSE
)
uniprot_small$uniprot_has_signal_peptide <- nzchar(trimws(uniprot_small$uniprot_signal_peptide))
uniprot_small$uniprot_explicitly_secreted <- grepl("SUBCELLULAR LOCATION:.*Secreted", uniprot_small$uniprot_subcellular_location, ignore.case = TRUE)
uniprot_small$uniprot_extracellular_go <- grepl("extracellular", uniprot_small$uniprot_go_cellular_component, ignore.case = TRUE)
uniprot_small$uniprot_platelet_or_blood_annotation <-
  grepl("platelet|erythrocyte|leukocyte|lymphocyte|neutrophil|monocyte|megakaryocyte", uniprot_small$uniprot_tissue_specificity, ignore.case = TRUE) |
  grepl("platelet alpha granule", paste(uniprot_small$uniprot_subcellular_location, uniprot_small$uniprot_go_cellular_component), ignore.case = TRUE)
uniprot_small$uniprot_liver_annotation <- grepl("liver|hepatocyte|hepatic", uniprot_small$uniprot_tissue_specificity, ignore.case = TRUE)

annotated <- merge(plasma, tissue_summary, by = "ensembl_id", all.x = TRUE, sort = FALSE)
annotated <- merge(annotated, immune_summary, by = "ensembl_id", all.x = TRUE, sort = FALSE)
annotated <- merge(annotated, uniprot_small, by = "accession", all.x = TRUE, sort = FALSE)
annotated <- annotated[match(plasma$accession, annotated$accession), , drop = FALSE]

ratio <- function(numerator, denominator) numerator / pmax(denominator, 0.1)
annotated$hpa_liver_to_skeletal_ratio <- ratio(annotated$hpa_liver_ntpm, annotated$hpa_skeletal_muscle_ntpm)
annotated$hpa_immune_to_skeletal_ratio <- ratio(annotated$hpa_immune_max_ntpm, annotated$hpa_skeletal_muscle_ntpm)
systemic_priority <- function(x, liver_min, liver_ratio, liver_high,
                              lymphoid_min, lymphoid_ratio,
                              immune_min, immune_ratio, immune_high, immune_high_ratio) {
  liver <- with(x,
    (hpa_tissue_max == "liver" & hpa_liver_ntpm >= liver_min &
      hpa_liver_ntpm >= liver_ratio * pmax(hpa_tissue_second_ntpm, 0.1)) |
      (hpa_liver_ntpm >= liver_high & hpa_liver_to_skeletal_ratio >= liver_ratio)
  )
  lymphoid <- with(x,
    hpa_tissue_max %in% c("bone marrow", "spleen", "lymph node", "tonsil", "thymus", "appendix") &
      hpa_tissue_max_ntpm >= lymphoid_min &
      hpa_tissue_max_ntpm >= lymphoid_ratio * pmax(hpa_skeletal_muscle_ntpm, 0.1)
  )
  immune <- with(x,
    (hpa_immune_max_ntpm >= immune_min & hpa_immune_to_skeletal_ratio >= immune_ratio) |
      (hpa_immune_max_ntpm >= immune_high & hpa_immune_to_skeletal_ratio >= immune_high_ratio)
  )
  list(liver = liver, lymphoid = lymphoid, immune = immune)
}

primary_priority <- systemic_priority(
  annotated, liver_min = 10, liver_ratio = 4, liver_high = 100,
  lymphoid_min = 10, lymphoid_ratio = 4,
  immune_min = 100, immune_ratio = 20, immune_high = 1000, immune_high_ratio = 5
)
annotated$hpa_liver_priority <- primary_priority$liver
annotated$hpa_lymphoid_priority <- primary_priority$lymphoid
annotated$hpa_immune_priority <- primary_priority$immune
annotated$competing_systemic_source <- with(annotated,
  hpa_liver_priority | hpa_lymphoid_priority | hpa_immune_priority |
    uniprot_platelet_or_blood_annotation | uniprot_liver_annotation
)
annotated$paraspinal_localized <- annotated$paraspinal_expression_support %in% c("strong_celltype_localization", "moderate_celltype_localization")
annotated$classical_or_explicit_secretory_support <- annotated$uniprot_has_signal_peptide | annotated$uniprot_explicitly_secreted

assign_origin_class <- function(x, systemic_source) {
  out <- rep("4_broad_or_unresolved", nrow(x))
  low_or_absent <- x$paraspinal_expression_support %in%
    c("not_detected_in_sc_reference", "low_or_sparse_expression") |
    is.na(x$paraspinal_expression_support)
  out[low_or_absent] <- "3_plasma_detected_paraspinal_low_or_absent"
  out[x$paraspinal_localized & !x$classical_or_explicit_secretory_support] <-
    "2_paraspinal_expressed_intracellular_or_damage_release"
  out[x$paraspinal_localized & x$classical_or_explicit_secretory_support & !systemic_source] <-
    "1_paraspinal_secretory_compatible_not_proven"
  out[systemic_source] <- "5_liver_or_hematologic_source_priority"
  out
}
annotated$cell_origin_evidence_class <- assign_origin_class(annotated, annotated$competing_systemic_source)

annotated$origin_interpretation <- c(
  "1_paraspinal_secretory_compatible_not_proven" = "Paraspinal cell-type localization plus secretion/extracellular annotation; compatible with, but does not prove, tissue contribution to plasma",
  "2_paraspinal_expressed_intracellular_or_damage_release" = "Paraspinal expression without classical secretion/extracellular support; intracellular or damage-release interpretation",
  "3_plasma_detected_paraspinal_low_or_absent" = "Detected in plasma export but little or no expression in the paraspinal reference",
  "4_broad_or_unresolved" = "Broad or insufficient localization; cellular source unresolved",
  "5_liver_or_hematologic_source_priority" = "HPA/UniProt evidence prioritizes liver or hematologic origin over a paraspinal source"
)[annotated$cell_origin_evidence_class]

write.table(annotated, file.path(out_dir, "plasma_protein_origin_evidence_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

class_counts <- as.data.frame(table(annotated$cell_origin_evidence_class), stringsAsFactors = FALSE)
names(class_counts) <- c("cell_origin_evidence_class", "n_proteins")
class_counts$percent <- 100 * class_counts$n_proteins / nrow(annotated)
write.table(class_counts, file.path(out_dir, "origin_evidence_class_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# HPA cut-offs are heuristic source-prioritization rules. Bracket the primary
# rules with permissive and stringent scenarios and report their consequences.
# Curated UniProt liver/blood annotations remain fixed across scenarios.
threshold_scenarios <- list(
  permissive = c(5, 3, 50, 5, 3, 50, 10, 500, 3),
  primary = c(10, 4, 100, 10, 4, 100, 20, 1000, 5),
  stringent = c(20, 6, 200, 20, 6, 200, 30, 2000, 8)
)
class_levels <- c(
  "1_paraspinal_secretory_compatible_not_proven",
  "2_paraspinal_expressed_intracellular_or_damage_release",
  "3_plasma_detected_paraspinal_low_or_absent",
  "4_broad_or_unresolved",
  "5_liver_or_hematologic_source_priority"
)
threshold_results <- do.call(rbind, lapply(names(threshold_scenarios), function(scenario) {
  z <- threshold_scenarios[[scenario]]
  p <- systemic_priority(
    annotated, liver_min = z[1], liver_ratio = z[2], liver_high = z[3],
    lymphoid_min = z[4], lymphoid_ratio = z[5], immune_min = z[6],
    immune_ratio = z[7], immune_high = z[8], immune_high_ratio = z[9]
  )
  systemic <- p$liver | p$lymphoid | p$immune |
    annotated$uniprot_platelet_or_blood_annotation | annotated$uniprot_liver_annotation
  cls <- assign_origin_class(annotated, systemic)
  counts <- table(factor(cls, levels = class_levels))
  data.frame(
    scenario = scenario,
    cell_origin_evidence_class = names(counts),
    n_proteins = as.integer(counts),
    n_high_confidence_class1 = sum(cls == class_levels[1] & annotated$fdr_confidence == "High"),
    stringsAsFactors = FALSE
  )
}))
write.table(
  threshold_results,
  file.path(out_dir, "source_priority_threshold_sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

candidate <- annotated[
  annotated$cell_origin_evidence_class == "1_paraspinal_secretory_compatible_not_proven" &
    annotated$fdr_confidence == "High",
  , drop = FALSE
]
candidate <- candidate[order(-candidate$dominant_cpm, -candidate$pooled_prevalence_in_dominant_type), , drop = FALSE]
write.table(candidate, file.path(out_dir, "high_confidence_paraspinal_secretory_compatible.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

pdf(file.path(out_dir, "origin_evidence_class_counts.pdf"), width = 9, height = 5.5)
par(mar = c(10, 5, 2, 1))
barplot(
  class_counts$n_proteins,
  names.arg = sub("^[1-5]_", "", class_counts$cell_origin_evidence_class),
  las = 2, ylab = "Proteins", col = c("#277DA1", "#577590", "#F9C74F", "#90BE6D", "#F94144"),
  main = "Origin-evidence classes (not causal source assignments)"
)
dev.off()

png(file.path(out_dir, "origin_evidence_class_counts.png"), width = 1800, height = 1100, res = 200)
par(mar = c(10, 5, 2, 1))
barplot(
  class_counts$n_proteins,
  names.arg = sub("^[1-5]_", "", class_counts$cell_origin_evidence_class),
  las = 2, ylab = "Proteins", col = c("#277DA1", "#577590", "#F9C74F", "#90BE6D", "#F94144"),
  main = "Origin-evidence classes (not causal source assignments)"
)
dev.off()

manifest <- data.frame(
  artifact = c("plasma_protein_origin_evidence_matrix.tsv", "origin_evidence_class_counts.tsv", "high_confidence_paraspinal_secretory_compatible.tsv"),
  n_rows = c(nrow(annotated), nrow(class_counts), nrow(candidate)),
  claim_boundary = c(
    "Evidence integration only; no protein is assigned a proven tissue source",
    "Descriptive evidence tiers, not source probabilities",
    "Candidates are compatible with paraspinal contribution but remain non-causal"
  )
)
write.table(manifest, file.path(out_dir, "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Integrated origin annotations for", nrow(annotated), "proteins\n")
print(class_counts)
cat("High-confidence paraspinal-secretory-compatible candidates:", nrow(candidate), "\n")
