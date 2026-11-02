#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(edgeR)
  library(limma)
  library(msigdbr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
input_xlsx <- if (length(args) >= 1) args[[1]] else "/Volumes/wang/paper14/ais_paraspinal_ecology/raw/bulk/GSE254300_readcount.xlsx"
bulk_dir <- if (length(args) >= 2) args[[2]] else "results/tissue/gse254300_bulk"
specificity_file <- if (length(args) >= 3) args[[3]] else "results/tissue/sc_reference/gene_celltype_specificity.tsv"
equivalence_file <- if (length(args) >= 4) args[[4]] else file.path(bulk_dir, "cross_etiology_equivalence", "gene_cross_etiology_equivalence.tsv.gz")
outdir <- if (length(args) >= 5) args[[5]] else "results/integration/composition_precision_sensitivity"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- as.data.frame(read_excel(input_xlsx, sheet = 1))
genes <- make.unique(as.character(dat$GeneID))
counts <- as.matrix(dat[, -1, drop = FALSE])
storage.mode(counts) <- "integer"
rownames(counts) <- genes
sample_id <- colnames(counts)
meta <- data.frame(
  sample_id = sample_id,
  group = sub("([A-Z]+).*", "\\1", sample_id),
  donor_number = sub("^[A-Z]+([0-9]+)[AT]$", "\\1", sample_id),
  side_code = sub(".*([AT])$", "\\1", sample_id),
  stringsAsFactors = FALSE
)
meta$side <- ifelse(meta$side_code == "A", "concave", "convex")
meta$donor <- paste(meta$group, meta$donor_number, sep = "_")

y <- DGEList(counts = counts)
keep <- filterByExpr(y, design = model.matrix(~0 + factor(meta$group):factor(meta$side)), min.count = 10)
y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
donor_design <- model.matrix(~0 + factor(meta$donor))
design <- cbind(
  donor_design,
  sideT = as.integer(meta$side == "convex"),
  interaction_IS_T = as.integer(meta$group == "IS" & meta$side == "convex")
)
contrast_shared <- c(rep(0, ncol(donor_design)), -1, -0.5)
y <- estimateDisp(y, design, robust = TRUE)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
sets <- split(
  c(hallmark$gene_symbol, reactome$gene_symbol),
  c(paste0("HALLMARK__", hallmark$gs_name), paste0("REACTOME__", reactome$gs_name))
)
sets <- lapply(sets, unique)

specificity <- read.delim(specificity_file, check.names = FALSE)
specificity <- specificity[is.finite(specificity$specificity_fraction) & specificity$dominant_cpm > 0, ]
specificity <- specificity[order(specificity$gene_symbol, -specificity$dominant_cpm), ]
specificity <- specificity[!duplicated(specificity$gene_symbol), ]
specificity_fraction <- setNames(specificity$specificity_fraction, specificity$gene_symbol)

run_camera <- function(label, eligible_genes) {
  # CAMERA is competitive: excluded genes must be removed from both the gene
  # sets and the background expression universe.
  eligible_genes <- intersect(rownames(y), eligible_genes)
  # Preserve the original TMM offsets so this sensitivity changes only the
  # competitive gene universe rather than also changing library scaling.
  y_filtered <- y[eligible_genes, , keep.lib.sizes = TRUE]
  index <- ids2indices(sets, rownames(y_filtered), remove.empty = TRUE)
  result <- edgeR::camera.DGEList(
    y_filtered, index = index, design = design, contrast = contrast_shared,
    use.ranks = TRUE, inter.gene.cor = 0.01
  )
  result$pathway <- rownames(result)
  result$gene_filter <- label
  result[, c("gene_filter", "pathway", "NGenes", "Direction", "PValue", "FDR")]
}

bulk_genes <- rownames(y)
broad_0p70 <- bulk_genes[is.na(specificity_fraction[bulk_genes]) | specificity_fraction[bulk_genes] < 0.70]
broad_0p50 <- bulk_genes[is.na(specificity_fraction[bulk_genes]) | specificity_fraction[bulk_genes] < 0.50]
camera_all <- run_camera("all_retained_genes", bulk_genes)
camera_0p70 <- run_camera("exclude_specificity_fraction_ge_0p70", broad_0p70)
camera_0p50 <- run_camera("exclude_specificity_fraction_ge_0p50", broad_0p50)
camera_sensitivity <- rbind(camera_all, camera_0p70, camera_0p50)
write.table(camera_sensitivity, file.path(outdir, "camera_celltype_specificity_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

primary_significant <- camera_all$pathway[camera_all$FDR < 0.05]
primary_direction <- setNames(camera_all$Direction, camera_all$pathway)
camera_summary <- do.call(rbind, lapply(split(camera_sensitivity, camera_sensitivity$gene_filter), function(x) {
  target <- x[x$pathway %in% primary_significant, ]
  data.frame(
    gene_filter = x$gene_filter[1],
    eligible_bulk_genes = switch(
      x$gene_filter[1],
      all_retained_genes = length(bulk_genes),
      exclude_specificity_fraction_ge_0p70 = length(broad_0p70),
      exclude_specificity_fraction_ge_0p50 = length(broad_0p50)
    ),
    original_significant_pathways = length(primary_significant),
    retained_fdr05 = sum(target$FDR < 0.05),
    retained_fdr10 = sum(target$FDR < 0.10),
    same_direction = sum(target$Direction == primary_direction[target$pathway]),
    stringsAsFactors = FALSE
  )
}))
write.table(camera_summary, file.path(outdir, "camera_celltype_specificity_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Test whether pathways are enriched among strictly equivalent genes using
# stratified label permutations that preserve abundance, effect-size and
# interaction-precision distributions.
equiv <- read.delim(gzfile(equivalence_file), check.names = FALSE)
shared <- equiv[equiv$shared_fdr05, ]
shared$strict <- as.integer(shared$strictly_conserved_margin_0p50)
shared$abs_shared_logFC <- abs(shared$shared_logFC)
shared$log_interaction_se <- log(pmax(shared$interaction_standard_error, 1e-8))
rank_bin <- function(x, bins = 3L) {
  pmin(bins, pmax(1L, ceiling(rank(x, ties.method = "average") / length(x) * bins)))
}
shared$stratum <- interaction(
  rank_bin(shared$shared_logCPM),
  rank_bin(shared$abs_shared_logFC),
  rank_bin(shared$log_interaction_se),
  drop = TRUE
)

membership_list <- lapply(sets, function(g) shared$gene_symbol %in% g)
membership_list <- membership_list[vapply(membership_list, function(x) sum(x) >= 10 && sum(x) <= nrow(shared) - 10, logical(1))]
membership <- do.call(cbind, membership_list)
storage.mode(membership) <- "numeric"
observed <- as.numeric(crossprod(membership, shared$strict))

permutations <- 5000L
set.seed(14016)
permuted_labels <- matrix(0, nrow = nrow(shared), ncol = permutations)
for (indices in split(seq_len(nrow(shared)), shared$stratum)) {
  strict_n <- sum(shared$strict[indices])
  if (strict_n == 0L) next
  if (strict_n == length(indices)) {
    permuted_labels[indices, ] <- 1
  } else {
    for (b in seq_len(permutations)) {
      permuted_labels[sample(indices, strict_n), b] <- 1
    }
  }
}
null_overlap <- crossprod(membership, permuted_labels)
expected <- rowMeans(null_overlap)
p_enriched <- (1 + rowSums(null_overlap >= observed)) / (permutations + 1)
p_depleted <- (1 + rowSums(null_overlap <= observed)) / (permutations + 1)
p_two_sided <- pmin(1, 2 * pmin(p_enriched, p_depleted))

pathway_total <- colSums(membership)
strict_total <- sum(shared$strict)
not_strict_total <- nrow(shared) - strict_total
a <- observed
b <- pathway_total - a
c <- strict_total - a
d <- not_strict_total - b
precision <- data.frame(
  pathway = colnames(membership),
  shared_genes_in_pathway = pathway_total,
  strict_genes_in_pathway = observed,
  expected_strict_genes_matched = expected,
  matched_enrichment_ratio = observed / pmax(expected, 1e-8),
  unadjusted_odds_ratio_haldane = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)),
  empirical_enrichment_PValue = p_enriched,
  empirical_depletion_PValue = p_depleted,
  empirical_two_sided_PValue = p_two_sided,
  direction = ifelse(observed >= expected, "enriched", "depleted"),
  stringsAsFactors = FALSE
)
precision$FDR <- p.adjust(precision$empirical_two_sided_PValue, method = "BH")
precision <- precision[order(precision$FDR, -abs(log2(pmax(precision$matched_enrichment_ratio, 1e-8)))), ]
write.table(precision, file.path(outdir, "strict_equivalence_pathway_precision_adjusted.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_metrics <- data.frame(
  metric = c(
    "original_camera_pathways_fdr05",
    "retained_after_excluding_specificity_ge_0p70_fdr05",
    "retained_after_excluding_specificity_ge_0p50_fdr05",
    "precision_adjusted_pathways_tested",
    "precision_adjusted_pathways_fdr05",
    "precision_adjusted_enriched_pathways_fdr05",
    "precision_adjusted_depleted_pathways_fdr05"
  ),
  value = c(
    length(primary_significant),
    camera_summary$retained_fdr05[camera_summary$gene_filter == "exclude_specificity_fraction_ge_0p70"],
    camera_summary$retained_fdr05[camera_summary$gene_filter == "exclude_specificity_fraction_ge_0p50"],
    nrow(precision), sum(precision$FDR < 0.05, na.rm = TRUE),
    sum(precision$FDR < 0.05 & precision$direction == "enriched", na.rm = TRUE),
    sum(precision$FDR < 0.05 & precision$direction == "depleted", na.rm = TRUE)
  )
)
write.table(summary_metrics, file.path(outdir, "analysis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_dat <- camera_sensitivity[camera_sensitivity$pathway %in% primary_significant, ]
plot_dat$label <- sub("^(HALLMARK|REACTOME)__", "", plot_dat$pathway)
plot_dat$score <- -log10(pmax(plot_dat$FDR, 1e-300)) * ifelse(plot_dat$Direction == "Up", 1, -1)
top <- camera_all[camera_all$FDR < 0.05, ]
top <- head(top[order(top$FDR), "pathway"], 15)
plot_dat <- plot_dat[plot_dat$pathway %in% top, ]
plot_dat$gene_filter <- factor(
  plot_dat$gene_filter,
  levels = c("all_retained_genes", "exclude_specificity_fraction_ge_0p70", "exclude_specificity_fraction_ge_0p50"),
  labels = c("All genes", "Exclude specificity >=0.70", "Exclude specificity >=0.50")
)
p <- ggplot(plot_dat, aes(gene_filter, label, fill = score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "Signed\n-log10 FDR")
ggsave(file.path(outdir, "celltype_specificity_pathway_sensitivity.pdf"), p, width = 8.2, height = 6.5)
ggsave(file.path(outdir, "celltype_specificity_pathway_sensitivity.png"), p, width = 8.2, height = 6.5, dpi = 240)
write.table(
  plot_dat[, c("pathway", "label", "gene_filter", "NGenes", "Direction", "PValue", "FDR", "score")],
  file.path(outdir, "source_data_figure4a.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

writeLines(c(
  "question\twhether shared pathways persist after removing genes with high single-cell cell-type specificity",
  "specificity_metric\tfraction of pooled reference expression contributed by the dominant cell type",
  "thresholds\texclude specificity fraction >=0.70 and >=0.50; genes absent from the single-cell reference retained",
  "pathway_test\tCAMERA competitive rank test after subsetting expression universe and gene-set indices to eligible genes; inter-gene correlation 0.01",
  paste0("precision_control\t", permutations, " strict-label permutations stratified by tertiles of bulk abundance, absolute shared effect and interaction standard error within 743 shared hits"),
  "boundary\tmarker exclusion reduces but cannot eliminate cell-composition confounding; the single-cell reference has two AIS donors",
  paste0("MSigDB_version\t", paste(unique(c(hallmark$db_version, reactome$db_version)), collapse = ";"))
), file.path(outdir, "analysis_manifest.tsv"))

cat("Composition and precision sensitivity analysis complete\n")
print(summary_metrics)
print(head(precision, 10))
