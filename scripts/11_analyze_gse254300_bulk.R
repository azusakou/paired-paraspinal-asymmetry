#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
input_xlsx <- if (length(args) >= 1) args[[1]] else "/Volumes/wang/paper14/ais_paraspinal_ecology/raw/bulk/GSE254300_readcount.xlsx"
outdir <- if (length(args) >= 2) args[[2]] else "results/tissue/gse254300_bulk"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- as.data.frame(read_excel(input_xlsx, sheet = 1))
stopifnot(colnames(dat)[1] == "GeneID")
genes <- make.unique(as.character(dat$GeneID))
counts <- as.matrix(dat[, -1, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- genes

sample_id <- colnames(counts)
stopifnot(all(grepl("^(IS|CS)[1-5][AT]$", sample_id)))
meta <- data.frame(
  sample_id = sample_id,
  group = sub("([A-Z]+).*", "\\1", sample_id),
  donor_number = sub("^[A-Z]+([0-9]+)[AT]$", "\\1", sample_id),
  side_code = sub(".*([AT])$", "\\1", sample_id),
  stringsAsFactors = FALSE
)
meta$side <- ifelse(meta$side_code == "A", "concave", "convex")
meta$donor <- paste(meta$group, meta$donor_number, sep = "_")
stopifnot(all(table(meta$donor) == 2))

y <- DGEList(counts = counts, samples = meta, genes = data.frame(gene_symbol = genes))
design_filter <- model.matrix(~ 0 + factor(meta$group):factor(meta$side))
keep <- filterByExpr(y, design = design_filter, min.count = 10)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

donor_design <- model.matrix(~ 0 + factor(meta$donor))
colnames(donor_design) <- paste0("donor_", levels(factor(meta$donor)))
sideT <- as.integer(meta$side == "convex")
interaction_IS_T <- as.integer(meta$group == "IS" & meta$side == "convex")
design <- cbind(donor_design, sideT = sideT, interaction_IS_T = interaction_IS_T)
stopifnot(qr(design)$rank == ncol(design))

y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

tests <- list(
  shared_concave_minus_convex = c(rep(0, ncol(donor_design)), -1, -0.5),
  AIS_concave_minus_convex = c(rep(0, ncol(donor_design)), -1, -1),
  CS_concave_minus_convex = c(rep(0, ncol(donor_design)), -1, 0),
  side_effect_interaction_AIS_minus_CS = c(rep(0, ncol(donor_design)), 0, -1)
)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
gene_sets <- split(
  c(hallmark$gene_symbol, reactome$gene_symbol),
  c(paste0("HALLMARK__", hallmark$gs_name), paste0("REACTOME__", reactome$gs_name))
)
gene_sets <- lapply(gene_sets, unique)
gene_set_indices <- ids2indices(gene_sets, rownames(y), remove.empty = TRUE)
camera_results <- lapply(names(tests), function(nm) {
  camera <- edgeR::camera.DGEList(
    y, index = gene_set_indices, design = design, contrast = tests[[nm]],
    use.ranks = TRUE, inter.gene.cor = 0.01
  )
  camera$pathway <- rownames(camera)
  camera$contrast <- nm
  camera <- camera[, c("contrast", "pathway", "NGenes", "Direction", "PValue", "FDR")]
  write.table(
    camera, file.path(outdir, paste0("camera_", nm, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  camera
})
names(camera_results) <- names(tests)
camera_shared <- camera_results$shared_concave_minus_convex
camera_by_etiology <- merge(
  camera_results$AIS_concave_minus_convex,
  camera_results$CS_concave_minus_convex,
  by = "pathway", suffixes = c("_AIS", "_CS")
)
camera_by_etiology$same_direction <- camera_by_etiology$Direction_AIS == camera_by_etiology$Direction_CS
camera_by_etiology$joint_fdr05 <- camera_by_etiology$FDR_AIS < 0.05 & camera_by_etiology$FDR_CS < 0.05
camera_by_etiology$mean_side_fdr05 <- camera_by_etiology$pathway %in% camera_shared$pathway[camera_shared$FDR < 0.05]
write.table(
  camera_by_etiology,
  file.path(outdir, "camera_AIS_CS_pathway_concordance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
primary_pathways <- camera_by_etiology[camera_by_etiology$mean_side_fdr05, ]
camera_pathway_summary <- data.frame(
  metric = c(
    "mean_side_camera_fdr05", "AIS_camera_fdr05", "CS_camera_fdr05",
    "joint_AIS_CS_camera_fdr05", "joint_AIS_CS_camera_fdr05_same_direction",
    "mean_side_hits_nominal_p05_in_both", "mean_side_hits_same_direction_in_both"
  ),
  value = c(
    sum(camera_shared$FDR < 0.05),
    sum(camera_results$AIS_concave_minus_convex$FDR < 0.05),
    sum(camera_results$CS_concave_minus_convex$FDR < 0.05),
    sum(camera_by_etiology$joint_fdr05),
    sum(camera_by_etiology$joint_fdr05 & camera_by_etiology$same_direction),
    sum(primary_pathways$PValue_AIS < 0.05 & primary_pathways$PValue_CS < 0.05),
    sum(primary_pathways$same_direction)
  )
)
write.table(
  camera_pathway_summary,
  file.path(outdir, "camera_cross_etiology_pathway_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

result_summary <- vector("list", length(tests))
names(result_summary) <- names(tests)
for (nm in names(tests)) {
  qlf <- glmQLFTest(fit, contrast = tests[[nm]])
  tab <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tab$gene_symbol <- rownames(tab)
  raw_tab <- qlf$table
  df_total <- qlf$df.total
  if (length(df_total) == 1L) df_total <- rep(df_total, nrow(raw_tab))
  standard_error <- abs(raw_tab$logFC) / sqrt(pmax(raw_tab$F, 0))
  standard_error[!is.finite(standard_error) | raw_tab$F <= 0] <- NA_real_
  uncertainty <- data.frame(
    gene_symbol = rownames(raw_tab),
    standard_error = standard_error,
    df_total = df_total,
    ci90_low = raw_tab$logFC - qt(0.95, df_total) * standard_error,
    ci90_high = raw_tab$logFC + qt(0.95, df_total) * standard_error,
    ci95_low = raw_tab$logFC - qt(0.975, df_total) * standard_error,
    ci95_high = raw_tab$logFC + qt(0.975, df_total) * standard_error,
    stringsAsFactors = FALSE
  )
  uncertainty <- uncertainty[match(tab$gene_symbol, uncertainty$gene_symbol), ]
  tab <- cbind(tab, uncertainty[, -1, drop = FALSE])
  tab$contrast <- nm
  tab$direction <- ifelse(tab$logFC > 0, "higher_concave_or_interaction", "lower_concave_or_interaction")
  tab <- tab[, c(
    "gene_symbol", "contrast", "logFC", "standard_error", "df_total",
    "ci90_low", "ci90_high", "ci95_low", "ci95_high",
    "logCPM", "F", "PValue", "FDR", "direction"
  )]
  write.table(tab, file.path(outdir, paste0(nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  result_summary[[nm]] <- data.frame(
    contrast = nm,
    tested_genes = nrow(tab),
    fdr_0.05 = sum(tab$FDR < 0.05),
    fdr_0.10 = sum(tab$FDR < 0.10),
    nominal_p_0.05 = sum(tab$PValue < 0.05),
    stringsAsFactors = FALSE
  )

  volcano <- ggplot(tab, aes(x = logFC, y = -log10(pmax(PValue, .Machine$double.xmin)), color = FDR < 0.05)) +
    geom_point(alpha = 0.55, size = 0.8) +
    scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#B2182B"), guide = "none") +
    geom_vline(xintercept = 0, linewidth = 0.25) +
    labs(title = gsub("_", " ", nm), x = "log2 fold change", y = "-log10(P)") +
    theme_bw(base_size = 9)
  ggsave(file.path(outdir, paste0(nm, "_volcano.pdf")), volcano, width = 4.5, height = 3.7)
}

summary_table <- do.call(rbind, result_summary)
write.table(summary_table, file.path(outdir, "contrast_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Leave one biological donor pair out at a time. Technical sides from the same
# donor are always removed together; normalization and dispersion are refit.
full_shared <- read.delim(file.path(outdir, "shared_concave_minus_convex.tsv"), check.names = FALSE)
loo_effects <- data.frame(gene_symbol = full_shared$gene_symbol, full_logFC = full_shared$logFC)
loo_correlations <- vector("list", length(unique(meta$donor)))
for (i in seq_along(unique(meta$donor))) {
  omitted <- unique(meta$donor)[i]
  use <- meta$donor != omitted
  meta_loo <- droplevels(meta[use, , drop = FALSE])
  y_loo <- y[, use, keep.lib.sizes = FALSE]
  y_loo$samples$norm.factors <- 1
  y_loo <- calcNormFactors(y_loo, method = "TMM")
  donor_loo <- model.matrix(~ 0 + factor(meta_loo$donor))
  design_loo <- cbind(
    donor_loo,
    sideT = as.integer(meta_loo$side == "convex"),
    interaction_IS_T = as.integer(meta_loo$group == "IS" & meta_loo$side == "convex")
  )
  stopifnot(qr(design_loo)$rank == ncol(design_loo))
  contrast_loo <- c(rep(0, ncol(donor_loo)), -1, -0.5)
  y_loo <- estimateDisp(y_loo, design_loo, robust = TRUE)
  fit_loo <- glmQLFit(y_loo, design_loo, robust = TRUE)
  tab_loo <- topTags(glmQLFTest(fit_loo, contrast = contrast_loo), n = Inf, sort.by = "none")$table
  lfc_loo <- tab_loo$logFC[match(full_shared$gene_symbol, rownames(tab_loo))]
  loo_effects[[paste0("omit_", omitted, "_logFC")]] <- lfc_loo
  loo_correlations[[i]] <- data.frame(
    omitted_donor = omitted,
    spearman_with_full = cor(full_shared$logFC, lfc_loo, method = "spearman", use = "complete.obs"),
    stringsAsFactors = FALSE
  )
}
loo_cols <- grep("^omit_.*_logFC$", names(loo_effects), value = TRUE)
loo_effects$all_loo_same_direction <- vapply(seq_len(nrow(loo_effects)), function(i) {
  all(sign(as.numeric(loo_effects[i, loo_cols])) == sign(loo_effects$full_logFC[i]))
}, logical(1))
loo_effects$loo_min_logFC <- apply(loo_effects[, loo_cols, drop = FALSE], 1, min, na.rm = TRUE)
loo_effects$loo_max_logFC <- apply(loo_effects[, loo_cols, drop = FALSE], 1, max, na.rm = TRUE)
write.table(loo_effects, gzfile(file.path(outdir, "leave_one_donor_out_shared_effects.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, loo_correlations), file.path(outdir, "leave_one_donor_out_correlations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
shared_idx <- full_shared$FDR < 0.05
write.table(
  data.frame(
    metric = c("full_shared_fdr05", "full_shared_fdr05_all_loo_same_direction", "minimum_loo_spearman", "median_loo_spearman"),
    value = c(
      sum(shared_idx), sum(shared_idx & loo_effects$all_loo_same_direction),
      min(vapply(loo_correlations, `[[`, numeric(1), "spearman_with_full")),
      median(vapply(loo_correlations, `[[`, numeric(1), "spearman_with_full"))
    )
  ),
  file.path(outdir, "leave_one_donor_out_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(meta, file.path(outdir, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  data.frame(gene_symbol = rownames(y), retained = TRUE, ave_log_cpm = aveLogCPM(y)),
  file.path(outdir, "retained_gene_universe.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

mds <- plotMDS(y, plot = FALSE)
mds_df <- data.frame(sample_id = sample_id, x = mds$x, y = mds$y, group = meta$group, side = meta$side)
write.table(mds_df, file.path(outdir, "mds_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p_mds <- ggplot(mds_df, aes(x = x, y = y, color = group, shape = side, label = sample_id)) +
  geom_point(size = 2.2) +
  geom_text(nudge_y = 0.08, size = 2.4, show.legend = FALSE) +
  theme_bw(base_size = 9) + labs(x = "Leading logFC dimension 1", y = "Leading logFC dimension 2")
ggsave(file.path(outdir, "mds.pdf"), p_mds, width = 5.2, height = 4.2)

manifest <- c(
  paste0("input_xlsx\t", normalizePath(input_xlsx)),
  paste0("input_genes\t", nrow(counts)),
  paste0("retained_genes\t", nrow(y)),
  "side_mapping\tA=concave;T=convex (GEO series annotation)",
  "model\tdonor fixed effects + convex side + AIS-by-convex interaction",
  "shared_side_contrast\tmean of AIS and congenital-scoliosis concave-minus-convex effects",
  "normalization\tTMM",
  "inference\tedgeR robust quasi-likelihood",
  paste0("R_version\t", R.version.string),
  paste0("edgeR_version\t", as.character(packageVersion("edgeR"))),
  paste0("limma_version\t", as.character(packageVersion("limma"))),
  paste0("MSigDB_version\t", paste(unique(c(hallmark$db_version, reactome$db_version)), collapse = ";"))
)
writeLines(manifest, file.path(outdir, "analysis_manifest.tsv"))
