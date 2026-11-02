#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
protein_features_file <- if (length(args) >= 1) args[[1]] else "results/audit/pxd024246/protein_features.tsv"
sc_specificity_file <- if (length(args) >= 2) args[[2]] else "results/tissue/sc_reference/gene_celltype_specificity.tsv"
sc_side_file <- if (length(args) >= 3) args[[3]] else "results/tissue/sc_reference/gene_celltype_side_effects.tsv"
bulk_dir <- if (length(args) >= 4) args[[4]] else "results/tissue/gse254300_bulk"
outdir <- if (length(args) >= 5) args[[5]] else "results/integration/tissue_support"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

proteins <- read.delim(protein_features_file, check.names = FALSE, stringsAsFactors = FALSE)
sc <- read.delim(sc_specificity_file, check.names = FALSE, stringsAsFactors = FALSE)
sc_side <- read.delim(sc_side_file, check.names = FALSE, stringsAsFactors = FALSE)
bulk_ais <- read.delim(file.path(bulk_dir, "AIS_concave_minus_convex.tsv"), stringsAsFactors = FALSE)
bulk_int <- read.delim(file.path(bulk_dir, "side_effect_interaction_AIS_minus_CS.tsv"), stringsAsFactors = FALSE)

extract_gene <- function(description) {
  hit <- regmatches(description, regexec("(?:^| )GN=([^ ]+)", description, perl = TRUE))
  vapply(hit, function(z) if (length(z) >= 2) z[[2]] else NA_character_, character(1))
}
proteins$gene_symbol <- extract_gene(proteins$description)
proteins$protein_order <- seq_len(nrow(proteins))

sc <- sc[!is.na(sc$gene_symbol) & nzchar(sc$gene_symbol), ]
sc <- sc[order(sc$gene_symbol, -sc$dominant_cpm, -sc$pooled_prevalence_in_dominant_type), ]
sc <- sc[!duplicated(sc$gene_symbol), ]

sc_keep <- sc[, c(
  "gene_symbol", "ensembl_id", "dominant_cell_type", "dominant_cpm",
  "specificity_fraction", "tau", "pooled_prevalence_in_dominant_type",
  "dominant_type_consistent_across_donors"
)]

integrated <- merge(proteins, sc_keep, by = "gene_symbol", all.x = TRUE, sort = FALSE)
integrated <- integrated[order(integrated$protein_order), ]

side_key <- paste(sc_side$gene_symbol, sc_side$cell_type, sep = "__")
target_key <- paste(integrated$gene_symbol, integrated$dominant_cell_type, sep = "__")
side_idx <- match(target_key, side_key)
integrated$sc_dominant_type_mean_log2cpm_concave_minus_convex <- sc_side$mean_log2cpm_concave_minus_convex[side_idx]
integrated$sc_dominant_type_side_direction_consistent <- sc_side$direction_consistent[side_idx]
integrated$sc_dominant_type_donors_positive <- sc_side$donors_positive[side_idx]

ais_idx <- match(integrated$gene_symbol, bulk_ais$gene_symbol)
int_idx <- match(integrated$gene_symbol, bulk_int$gene_symbol)
integrated$bulk_AIS_concave_minus_convex_logFC <- bulk_ais$logFC[ais_idx]
integrated$bulk_AIS_PValue <- bulk_ais$PValue[ais_idx]
integrated$bulk_AIS_FDR <- bulk_ais$FDR[ais_idx]
integrated$bulk_interaction_AIS_minus_CS_logFC <- bulk_int$logFC[int_idx]
integrated$bulk_interaction_PValue <- bulk_int$PValue[int_idx]
integrated$bulk_interaction_FDR <- bulk_int$FDR[int_idx]

has_sc <- !is.na(integrated$dominant_cpm)
strong <- has_sc & integrated$dominant_cpm >= 5 &
  integrated$pooled_prevalence_in_dominant_type >= 0.10 &
  integrated$specificity_fraction >= 0.50 &
  integrated$dominant_type_consistent_across_donors
moderate <- has_sc & integrated$dominant_cpm >= 1 &
  integrated$pooled_prevalence_in_dominant_type >= 0.02 &
  integrated$specificity_fraction >= 0.40 &
  integrated$dominant_type_consistent_across_donors
broad <- has_sc & integrated$dominant_cpm >= 1 &
  integrated$pooled_prevalence_in_dominant_type >= 0.02

integrated$paraspinal_expression_support <- "low_or_sparse_expression"
integrated$paraspinal_expression_support[!has_sc & !is.na(integrated$gene_symbol)] <- "not_detected_in_sc_reference"
integrated$paraspinal_expression_support[is.na(integrated$gene_symbol)] <- "unmapped_protein_export"
integrated$paraspinal_expression_support[broad] <- "broad_expression"
integrated$paraspinal_expression_support[moderate] <- "moderate_celltype_localization"
integrated$paraspinal_expression_support[strong] <- "strong_celltype_localization"

integrated$bulk_AIS_side_support <- !is.na(integrated$bulk_AIS_FDR) & integrated$bulk_AIS_FDR < 0.05
integrated$bulk_AIS_side_suggestive <- !is.na(integrated$bulk_AIS_PValue) &
  integrated$bulk_AIS_PValue < 0.05 & abs(integrated$bulk_AIS_concave_minus_convex_logFC) >= 0.5
integrated$sc_bulk_side_direction_concordant <- with(integrated,
  !is.na(sc_dominant_type_mean_log2cpm_concave_minus_convex) &
    !is.na(bulk_AIS_concave_minus_convex_logFC) &
    sign(sc_dominant_type_mean_log2cpm_concave_minus_convex) == sign(bulk_AIS_concave_minus_convex_logFC)
)

integrated$origin_claim_status <- "expression_only_not_origin"
integrated$origin_claim_status[integrated$paraspinal_expression_support == "not_detected_in_sc_reference"] <- "no_sc_expression_support"
integrated$origin_claim_status[integrated$paraspinal_expression_support == "unmapped_protein_export"] <- "unmapped"

write.table(integrated, file.path(outdir, "plasma_protein_tissue_support_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_tab <- as.data.frame(table(
  fdr_confidence = integrated$fdr_confidence,
  expression_support = integrated$paraspinal_expression_support,
  useNA = "ifany"
))
summary_tab <- summary_tab[summary_tab$Freq > 0, ]
write.table(summary_tab, file.path(outdir, "support_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_dat <- integrated[
  integrated$fdr_confidence == "High" & has_sc &
    integrated$dominant_cpm >= 1 & integrated$pooled_prevalence_in_dominant_type >= 0.02,
]
p <- ggplot(plot_dat, aes(
  x = specificity_fraction,
  y = pooled_prevalence_in_dominant_type,
  color = dominant_cell_type,
  size = log10(dominant_cpm + 1)
)) +
  geom_point(alpha = 0.72) +
  geom_vline(xintercept = c(0.4, 0.5), linetype = c(3, 2), linewidth = 0.3) +
  geom_hline(yintercept = c(0.02, 0.10), linetype = c(3, 2), linewidth = 0.3) +
  labs(
    x = "Cell-type specificity fraction",
    y = "Expression prevalence in dominant cell type",
    color = "Dominant cell type",
    size = "log10(CPM + 1)"
  ) +
  theme_bw(base_size = 9)
ggsave(file.path(outdir, "protein_gene_sc_expression_support.pdf"), p, width = 6.2, height = 4.5)

manifest <- c(
  paste0("plasma_proteins\t", nrow(integrated)),
  paste0("proteins_with_embedded_gene_symbol\t", sum(!is.na(integrated$gene_symbol))),
  paste0("proteins_with_sc_match\t", sum(has_sc)),
  paste0("strong_celltype_localization\t", sum(strong)),
  paste0("moderate_celltype_localization_including_strong\t", sum(moderate)),
  "strong_rule\tdominant CPM >=5; prevalence >=0.10; specificity fraction >=0.50; same dominant type in both donors",
  "moderate_rule\tdominant CPM >=1; prevalence >=0.02; specificity fraction >=0.40; same dominant type in both donors",
  "claim_boundary\tThese tiers quantify paraspinal expression localization only. They do not establish secretion or plasma cell of origin.",
  paste0("R_version\t", R.version.string)
)
writeLines(manifest, file.path(outdir, "analysis_manifest.tsv"))

