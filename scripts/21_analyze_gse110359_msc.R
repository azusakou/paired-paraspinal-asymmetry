#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
matrix_file <- if (length(args) >= 1) args[[1]] else "/Volumes/wang/paper15/ais_mesenchymal_axis/raw/gse110359/GSE110359_series_matrix.txt.gz"
annotation_file <- if (length(args) >= 2) args[[2]] else "metadata/gse110359/probe_annotation.tsv"
sample_file <- if (length(args) >= 3) args[[3]] else "metadata/gse110359/sample_metadata.tsv"
outdir <- if (length(args) >= 4) args[[4]] else "results/tissue/gse110359_msc"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

lines <- readLines(gzfile(matrix_file), warn = FALSE)
start <- which(lines == "!series_matrix_table_begin") + 1L
end <- which(lines == "!series_matrix_table_end") - 1L
stopifnot(length(start) == 1L, length(end) == 1L, end > start)
expr_df <- read.delim(textConnection(lines[start:end]), check.names = FALSE, quote = "\"")
rm(lines)
probe_ids <- expr_df[[1]]
expr <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- probe_ids

annotation <- read.delim(annotation_file, check.names = FALSE)
samples <- read.delim(sample_file, check.names = FALSE)
stopifnot(identical(colnames(expr), samples$sample_id))
idx <- match(rownames(expr), annotation$probe_id)
stopifnot(!anyNA(idx))
annotation <- annotation[idx, , drop = FALSE]

keep <- annotation$mapping_status == "unique" & nzchar(annotation$gene_symbol)
expr <- expr[keep, , drop = FALSE]
annotation <- annotation[keep, , drop = FALSE]

# Collapse probes without using disease labels: retain the highest-IQR probe per gene.
probe_iqr <- apply(expr, 1, IQR, na.rm = TRUE)
probe_order <- order(annotation$gene_symbol, -probe_iqr, rownames(expr))
representative <- probe_order[!duplicated(annotation$gene_symbol[probe_order])]
expr_gene <- expr[representative, , drop = FALSE]
gene_annotation <- annotation[representative, , drop = FALSE]
rownames(expr_gene) <- gene_annotation$gene_symbol

group <- factor(samples$disease_state, levels = c("healthy", "ais"))
stopifnot(all(table(group) == c(5, 12)))
design <- model.matrix(~ group)
colnames(design) <- c("intercept", "AIS_minus_healthy")

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
gene_sets <- split(
  c(hallmark$gene_symbol, reactome$gene_symbol),
  c(paste0("HALLMARK__", hallmark$gs_name), paste0("REACTOME__", reactome$gs_name))
)
gene_sets <- lapply(gene_sets, unique)
gene_set_indices <- ids2indices(gene_sets, rownames(expr_gene), remove.empty = TRUE)
camera_msc <- camera(
  expr_gene, index = gene_set_indices, design = design, contrast = 2,
  use.ranks = TRUE, inter.gene.cor = 0.01
)
camera_msc$pathway <- rownames(camera_msc)
camera_msc <- camera_msc[, c("pathway", "NGenes", "Direction", "PValue", "FDR")]
write.table(camera_msc, file.path(outdir, "camera_AIS_minus_healthy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

fit_one <- function(mat, method) {
  fit <- lmFit(mat, design)
  fit <- eBayes(fit, robust = TRUE, trend = TRUE)
  tab <- topTable(fit, coef = "AIS_minus_healthy", n = Inf, sort.by = "P")
  tab$gene_symbol <- rownames(tab)
  tab$normalization <- method
  tab$direction <- ifelse(tab$logFC > 0, "higher_AIS", "lower_AIS")
  tab[, c("gene_symbol", "normalization", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "direction")]
}

result_deposited <- fit_one(expr_gene, "deposited_RMA")
expr_quantile <- normalizeBetweenArrays(expr_gene, method = "quantile")
result_quantile <- fit_one(expr_quantile, "additional_quantile")
write.table(result_deposited, file.path(outdir, "AIS_minus_healthy_deposited_RMA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(result_quantile, file.path(outdir, "AIS_minus_healthy_additional_quantile.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sensitivity <- merge(
  result_deposited[, c("gene_symbol", "logFC", "adj.P.Val")],
  result_quantile[, c("gene_symbol", "logFC", "adj.P.Val")],
  by = "gene_symbol", suffixes = c("_deposited", "_quantile")
)
sensitivity$direction_concordant <- sign(sensitivity$logFC_deposited) == sign(sensitivity$logFC_quantile)
write.table(sensitivity, file.path(outdir, "normalization_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sample_qc <- data.frame(
  sample_id = samples$sample_id,
  disease_state = samples$disease_state,
  median_expression = apply(expr_gene, 2, median),
  expression_iqr = apply(expr_gene, 2, IQR),
  sample_correlation_median = apply(cor(expr_gene, method = "spearman"), 2, median),
  stringsAsFactors = FALSE
)
write.table(sample_qc, file.path(outdir, "sample_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pca <- prcomp(t(expr_gene), scale. = FALSE)
pca_df <- data.frame(
  sample_id = samples$sample_id,
  disease_state = samples$disease_state,
  PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)
write.table(pca_df, file.path(outdir, "pca_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p <- ggplot(pca_df, aes(PC1, PC2, color = disease_state, label = sample_id)) +
  geom_point(size = 2.4) + geom_text(nudge_y = 5, size = 2.3, show.legend = FALSE) +
  theme_bw(base_size = 9) + labs(title = "GSE110359 BM-MSC expression")
ggsave(file.path(outdir, "pca.pdf"), p, width = 5.2, height = 4.2)

summary <- data.frame(
  normalization = c("deposited_RMA", "additional_quantile"),
  tested_genes = c(nrow(result_deposited), nrow(result_quantile)),
  fdr_0.05 = c(sum(result_deposited$adj.P.Val < 0.05), sum(result_quantile$adj.P.Val < 0.05)),
  fdr_0.10 = c(sum(result_deposited$adj.P.Val < 0.10), sum(result_quantile$adj.P.Val < 0.10)),
  nominal_p_0.05 = c(sum(result_deposited$P.Value < 0.05), sum(result_quantile$P.Value < 0.05)),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(outdir, "analysis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  data.frame(
    metric = c("unique_mapped_probes", "representative_genes", "effect_spearman", "direction_concordance", "fdr05_overlap"),
    value = c(
      nrow(expr), nrow(expr_gene),
      cor(sensitivity$logFC_deposited, sensitivity$logFC_quantile, method = "spearman"),
      mean(sensitivity$direction_concordant),
      sum(sensitivity$adj.P.Val_deposited < 0.05 & sensitivity$adj.P.Val_quantile < 0.05)
    )
  ),
  file.path(outdir, "sensitivity_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- c(
  paste0("input_matrix\t", normalizePath(matrix_file)),
  paste0("input_annotation\t", normalizePath(annotation_file)),
  "contrast\tAIS minus healthy BM-MSCs",
  "probe_collapse\thighest IQR per uniquely mapped gene; disease labels not used",
  "primary_normalization\tdeposited RMA values",
  "sensitivity_normalization\tadditional quantile normalization",
  "inference\tlimma empirical Bayes, robust=TRUE, trend=TRUE",
  "covariates\tnone available at individual level; age only reported as <18 and sex unavailable",
  paste0("R_version\t", R.version.string),
  paste0("limma_version\t", as.character(packageVersion("limma")))
  ,paste0("MSigDB_version\t", paste(unique(c(hallmark$db_version, reactome$db_version)), collapse = ";"))
)
writeLines(manifest, file.path(outdir, "analysis_manifest.tsv"))

cat("GSE110359 gene-level analysis complete\n")
print(summary)
