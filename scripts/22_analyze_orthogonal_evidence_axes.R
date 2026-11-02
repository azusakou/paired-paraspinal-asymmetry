#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
bulk_dir <- if (length(args) >= 1) args[[1]] else "results/tissue/gse254300_bulk"
msc_dir <- if (length(args) >= 2) args[[2]] else "results/tissue/gse110359_msc"
sc_file <- if (length(args) >= 3) args[[3]] else "results/tissue/sc_reference/gene_celltype_side_effects.tsv"
origin_file <- if (length(args) >= 4) args[[4]] else "results/integration/origin_annotation/plasma_protein_origin_evidence_matrix.tsv"
outdir <- if (length(args) >= 5) args[[5]] else "results/integration/orthogonal_axes"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_contrast <- function(name) read.delim(file.path(bulk_dir, paste0(name, ".tsv")), check.names = FALSE)
shared <- read_contrast("shared_concave_minus_convex")
ais <- read_contrast("AIS_concave_minus_convex")
cs <- read_contrast("CS_concave_minus_convex")
interaction <- read_contrast("side_effect_interaction_AIS_minus_CS")
msc <- read.delim(file.path(msc_dir, "AIS_minus_healthy_deposited_RMA.tsv"), check.names = FALSE)
sc <- read.delim(sc_file, check.names = FALSE)
origin <- read.delim(origin_file, check.names = FALSE)

rename_cols <- function(x, prefix, cols) {
  keep <- x[, c("gene_symbol", cols), drop = FALSE]
  names(keep)[-1] <- paste0(prefix, "_", names(keep)[-1])
  keep
}

gene_map <- Reduce(function(x, y) merge(x, y, by = "gene_symbol", all = TRUE), list(
  rename_cols(shared, "shared", c("logFC", "logCPM", "F", "PValue", "FDR")),
  rename_cols(ais, "ais", c("logFC", "PValue", "FDR")),
  rename_cols(cs, "cs", c("logFC", "PValue", "FDR")),
  rename_cols(interaction, "interaction", c("logFC", "PValue", "FDR")),
  rename_cols(msc, "msc", c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"))
))

gene_map$shared_adaptation <- with(gene_map,
  shared_FDR < 0.05 & sign(ais_logFC) == sign(cs_logFC) & interaction_FDR >= 0.10
)
gene_map$ais_msc_fdr10 <- gene_map$msc_adj.P.Val < 0.10

sc_fap <- sc[sc$cell_type == "FAP", c(
  "gene_symbol", "mean_log2cpm_concave_minus_convex", "direction_consistent",
  "donors_positive", "donors_negative"
)]
names(sc_fap)[-1] <- paste0("sc_fap_", names(sc_fap)[-1])
gene_map <- merge(gene_map, sc_fap, by = "gene_symbol", all.x = TRUE)

origin$class_order <- match(origin$cell_origin_evidence_class, c(
  "1_paraspinal_secretory_compatible_not_proven",
  "2_paraspinal_expressed_intracellular_or_damage_release",
  "4_broad_or_unresolved",
  "3_plasma_detected_paraspinal_low_or_absent",
  "5_liver_or_hematologic_source_priority"
))
origin$confidence_order <- match(origin$fdr_confidence, c("High", "Medium", "Low"))
origin <- origin[order(origin$gene_symbol, origin$class_order, origin$confidence_order), ]
origin_gene <- origin[!is.na(origin$gene_symbol) & nzchar(origin$gene_symbol) & !duplicated(origin$gene_symbol), c(
  "gene_symbol", "accession", "fdr_confidence", "dominant_cell_type", "dominant_cpm",
  "specificity_fraction", "pooled_prevalence_in_dominant_type",
  "cell_origin_evidence_class", "classical_or_explicit_secretory_support",
  "competing_systemic_source"
)]
names(origin_gene)[-1] <- paste0("plasma_", names(origin_gene)[-1])
gene_map <- merge(gene_map, origin_gene, by = "gene_symbol", all.x = TRUE)
gene_map$plasma_detected <- !is.na(gene_map$plasma_accession)
gene_map$bridge_candidate <- with(gene_map,
  shared_adaptation & plasma_detected &
    plasma_cell_origin_evidence_class == "1_paraspinal_secretory_compatible_not_proven"
)
write.table(gene_map, file.path(outdir, "gene_evidence_axis_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

bridge <- gene_map[gene_map$bridge_candidate, , drop = FALSE]
bridge <- bridge[order(bridge$shared_FDR, -abs(bridge$shared_logFC)), ]
write.table(bridge, file.path(outdir, "shared_adaptation_plasma_bridge_candidates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

scatter <- gene_map[is.finite(gene_map$ais_logFC) & is.finite(gene_map$cs_logFC), ]
p_scatter <- ggplot(scatter, aes(ais_logFC, cs_logFC, color = shared_adaptation)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.35) +
  geom_point(alpha = 0.45, size = 0.8) +
  geom_point(data = scatter[scatter$bridge_candidate, ], shape = 21, size = 3, stroke = 0.8, fill = "#FEE08B") +
  geom_label_repel(
    data = scatter[scatter$bridge_candidate, ], aes(label = gene_symbol),
    size = 3, box.padding = 0.7, point.padding = 0.5, force = 4,
    fill = scales::alpha("white", 0.88), label.size = 0.15,
    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE,
    seed = 14016
  ) +
  scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "#2166AC"), na.value = "grey85", name = "Shared adaptation") +
  coord_equal() + theme_bw(base_size = 9) +
  labs(x = "AIS concave - convex log2 fold change", y = "Congenital scoliosis concave - convex log2 fold change")
ggsave(file.path(outdir, "AIS_vs_congenital_side_effects.pdf"), p_scatter, width = 6.3, height = 5.5)
ggsave(file.path(outdir, "AIS_vs_congenital_side_effects.png"), p_scatter, width = 6.3, height = 5.5, dpi = 240)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
msig_version <- unique(c(hallmark$db_version, reactome$db_version))
sets <- split(c(hallmark$gene_symbol, reactome$gene_symbol), c(paste0("HALLMARK__", hallmark$gs_name), paste0("REACTOME__", reactome$gs_name)))
sets <- lapply(sets, unique)

run_fgsea <- function(stats, axis) {
  stats <- stats[is.finite(stats) & !is.na(names(stats)) & nzchar(names(stats))]
  stats <- sort(stats[!duplicated(names(stats))], decreasing = TRUE)
  set.seed(14016)
  result <- as.data.frame(fgseaMultilevel(pathways = sets, stats = stats, minSize = 15, maxSize = 500, eps = 0))
  result$axis <- axis
  result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";", character(1))
  result[order(result$padj, -abs(result$NES)), c("axis", "pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge")]
}

bulk_stats <- setNames(sign(shared$logFC) * sqrt(pmax(shared$F, 0)), shared$gene_symbol)
msc_stats <- setNames(msc$t, msc$gene_symbol)
fap <- sc[sc$cell_type == "FAP" & sc$direction_consistent, ]
fap_stats <- setNames(fap$mean_log2cpm_concave_minus_convex, fap$gene_symbol)

fg_bulk <- run_fgsea(bulk_stats, "shared_side_adaptation")
fg_msc <- run_fgsea(msc_stats, "AIS_BM_MSC_disease")
fg_fap <- run_fgsea(fap_stats, "AIS_FAP_side_descriptive")
write.table(fg_bulk, file.path(outdir, "pathways_shared_side_adaptation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fg_msc, file.path(outdir, "pathways_AIS_BM_MSC_disease.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fg_fap, file.path(outdir, "pathways_AIS_FAP_side_descriptive.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pathway_map_fgsea <- merge(
  fg_bulk[, c("pathway", "NES", "padj")],
  fg_msc[, c("pathway", "NES", "padj")],
  by = "pathway", suffixes = c("_shared_side", "_msc")
)
pathway_map_fgsea$both_fdr05 <- pathway_map_fgsea$padj_shared_side < 0.05 & pathway_map_fgsea$padj_msc < 0.05
pathway_map_fgsea$direction_concordant <- sign(pathway_map_fgsea$NES_shared_side) == sign(pathway_map_fgsea$NES_msc)
pathway_map_fgsea$evidence_axis_class <- "not_jointly_significant"
pathway_map_fgsea$evidence_axis_class[pathway_map_fgsea$both_fdr05 & pathway_map_fgsea$direction_concordant] <- "convergent_disease_and_adaptation"
pathway_map_fgsea$evidence_axis_class[pathway_map_fgsea$both_fdr05 & !pathway_map_fgsea$direction_concordant] <- "opposed_disease_and_adaptation"
write.table(pathway_map_fgsea[order(pathway_map_fgsea$padj_shared_side, pathway_map_fgsea$padj_msc), ], file.path(outdir, "pathway_evidence_axis_map_fgsea_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# CAMERA is the confirmatory pathway analysis because it adjusts competitive
# gene-set tests for inter-gene correlation within each expression dataset.
camera_bulk <- read.delim(file.path(bulk_dir, "camera_shared_concave_minus_convex.tsv"), check.names = FALSE)
camera_msc <- read.delim(file.path(msc_dir, "camera_AIS_minus_healthy.tsv"), check.names = FALSE)
pathway_map <- merge(camera_bulk, camera_msc, by = "pathway", suffixes = c("_shared_side", "_msc"))
pathway_map$both_fdr05 <- pathway_map$FDR_shared_side < 0.05 & pathway_map$FDR_msc < 0.05
pathway_map$direction_concordant <- pathway_map$Direction_shared_side == pathway_map$Direction_msc
pathway_map$evidence_axis_class <- "not_jointly_significant"
pathway_map$evidence_axis_class[pathway_map$both_fdr05 & pathway_map$direction_concordant] <- "convergent_disease_and_adaptation"
pathway_map$evidence_axis_class[pathway_map$both_fdr05 & !pathway_map$direction_concordant] <- "opposed_disease_and_adaptation"
write.table(pathway_map[order(pathway_map$FDR_shared_side, pathway_map$FDR_msc), ], file.path(outdir, "pathway_evidence_axis_map_camera.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary <- data.frame(
  metric = c(
    "shared_adaptation_genes_fdr05", "shared_adaptation_same_direction",
    "interaction_genes_fdr10", "msc_genes_fdr05", "msc_genes_fdr10",
    "bridge_candidate_protein_entries", "bridge_candidate_unique_genes",
    "shared_pathways_camera_fdr05", "msc_pathways_camera_fdr05", "joint_concordant_pathways_camera_fdr05",
    "joint_opposed_pathways_camera_fdr05"
  ),
  value = c(
    sum(shared$FDR < 0.05),
    sum(shared$FDR < 0.05 & sign(ais$logFC[match(shared$gene_symbol, ais$gene_symbol)]) == sign(cs$logFC[match(shared$gene_symbol, cs$gene_symbol)])),
    sum(interaction$FDR < 0.10),
    sum(msc$adj.P.Val < 0.05), sum(msc$adj.P.Val < 0.10),
    nrow(bridge), length(unique(bridge$gene_symbol)),
    sum(camera_bulk$FDR < 0.05, na.rm = TRUE), sum(camera_msc$FDR < 0.05, na.rm = TRUE),
    sum(pathway_map$evidence_axis_class == "convergent_disease_and_adaptation"),
    sum(pathway_map$evidence_axis_class == "opposed_disease_and_adaptation")
  )
)
write.table(summary, file.path(outdir, "analysis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Confirmatory Figure 3 uses CAMERA, matching the manuscript's inferential
# pathway analysis. FGSEA remains a separately named sensitivity output.
camera_plot <- rbind(
  transform(camera_bulk[, c("pathway", "NGenes", "Direction", "PValue", "FDR")], axis = "Paraspinal mean side contrast"),
  transform(camera_msc[, c("pathway", "NGenes", "Direction", "PValue", "FDR")], axis = "AIS versus healthy BM-MSC")
)
camera_plot <- camera_plot[camera_plot$FDR < 0.05, ]
camera_plot <- camera_plot[order(camera_plot$FDR, camera_plot$PValue), ]
camera_plot <- do.call(rbind, lapply(split(camera_plot, camera_plot$axis), head, 10))
camera_plot$label <- gsub("_", " ", sub("^(HALLMARK|REACTOME)__", "", camera_plot$pathway))
camera_plot$label <- vapply(
  camera_plot$label,
  function(x) paste(strwrap(x, width = 52), collapse = "\n"),
  character(1)
)
camera_plot$signed_evidence <- -log10(pmax(camera_plot$FDR, 1e-300)) * ifelse(camera_plot$Direction == "Up", 1, -1)
camera_plot$axis <- factor(
  camera_plot$axis,
  levels = c("Paraspinal mean side contrast", "AIS versus healthy BM-MSC")
)
p_camera <- ggplot(camera_plot, aes(signed_evidence, reorder(label, signed_evidence))) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_point(size = 2.2, alpha = 0.9, color = "#2166AC") +
  facet_wrap(~axis, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 8.5) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 7.2, lineheight = 0.9),
    panel.grid.major.y = element_blank()
  ) +
  labs(x = "Signed -log10 CAMERA FDR", y = NULL)
ggsave(file.path(outdir, "orthogonal_axis_pathways_camera.pdf"), p_camera, width = 8.2, height = 8.6)
ggsave(file.path(outdir, "orthogonal_axis_pathways_camera.png"), p_camera, width = 8.2, height = 8.6, dpi = 300)
write.table(
  camera_plot[, c("axis", "pathway", "label", "NGenes", "Direction", "PValue", "FDR", "signed_evidence")],
  file.path(outdir, "source_data_figure3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

fgsea_plot <- rbind(fg_bulk, fg_msc)
fgsea_plot <- fgsea_plot[fgsea_plot$padj < 0.05, ]
fgsea_plot <- fgsea_plot[order(fgsea_plot$padj, -abs(fgsea_plot$NES)), ]
fgsea_plot <- do.call(rbind, lapply(split(fgsea_plot, fgsea_plot$axis), head, 12))
fgsea_plot$label <- sub("^(HALLMARK|REACTOME)__", "", fgsea_plot$pathway)
p_fgsea <- ggplot(fgsea_plot, aes(NES, reorder(label, NES), size = -log10(pmax(padj, 1e-300)), color = axis)) +
  geom_point(alpha = 0.85) + facet_wrap(~axis, scales = "free_y") +
  theme_bw(base_size = 9) + labs(x = "Normalized enrichment score", y = NULL, size = "-log10 FDR")
ggsave(file.path(outdir, "orthogonal_axis_pathways_fgsea_sensitivity.pdf"), p_fgsea, width = 10, height = 7)
ggsave(file.path(outdir, "orthogonal_axis_pathways_fgsea_sensitivity.png"), p_fgsea, width = 10, height = 7, dpi = 300)

manifest <- c(
  "concept\torthogonal contrasts separate shared curvature-associated adaptation from AIS-associated BM-MSC state",
  "shared_adaptation_definition\tbulk shared-side FDR<0.05; AIS and CS same direction; interaction FDR>=0.10",
  "bridge_definition\tshared adaptation plus PXD detection plus paraspinal-secretory-compatible source class",
  "causal_boundary\tevidence axes do not establish temporal or causal ordering",
  "confirmatory_pathway_test\tCAMERA competitive test with inter-gene correlation=0.01 and rank statistic",
  "pathway_sensitivity_test\tfgsea multilevel pre-ranked enrichment",
  paste0("MSigDB_version\t", paste(msig_version, collapse = ";")),
  paste0("fgsea_version\t", as.character(packageVersion("fgsea"))),
  paste0("msigdbr_version\t", as.character(packageVersion("msigdbr")))
)
writeLines(manifest, file.path(outdir, "analysis_manifest.tsv"))

cat("Orthogonal evidence-axis analysis complete\n")
print(summary)
