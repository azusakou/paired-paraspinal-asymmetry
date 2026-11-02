#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
matrix_tsv <- if (length(args) >= 1) args[[1]] else "results/audit/pxd024246/abundance_matrix_raw_export.tsv"
features_tsv <- if (length(args) >= 2) args[[2]] else "results/audit/pxd024246/protein_features.tsv"
outdir <- if (length(args) >= 3) args[[3]] else "results/proteomics/unsupervised_qc"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(matrix_tsv, check.names = FALSE, stringsAsFactors = FALSE)
features <- read.delim(features_tsv, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(dat$accession == features$accession)
x <- as.matrix(dat[, -1, drop = FALSE])
storage.mode(x) <- "numeric"
stopifnot(all(x[is.finite(x)] > 0))
log2x <- log2(x)

sample_metrics <- data.frame(
  abundance_column = colnames(x),
  n_observed = colSums(!is.na(x)),
  missing_fraction = colMeans(is.na(x)),
  median_log2_abundance = apply(log2x, 2, median, na.rm = TRUE),
  iqr_log2_abundance = apply(log2x, 2, IQR, na.rm = TRUE),
  stringsAsFactors = FALSE
)

protein_metrics <- data.frame(
  accession = dat$accession,
  fdr_confidence = features$fdr_confidence,
  n_observed = rowSums(!is.na(x)),
  missing_fraction = rowMeans(is.na(x)),
  median_log2_abundance = apply(log2x, 1, median, na.rm = TRUE),
  iqr_log2_abundance = apply(log2x, 1, IQR, na.rm = TRUE),
  stringsAsFactors = FALSE
)

complete <- rowSums(is.na(log2x)) == 0 & features$fdr_confidence == "High"
if (sum(complete) < 20) stop("Fewer than 20 high-confidence complete-case proteins; PCA aborted")
cc <- log2x[complete, , drop = FALSE]

median_centered <- sweep(cc, 2, apply(cc, 2, median), "-")
quantile_normalized <- normalizeBetweenArrays(cc, method = "quantile")
strategies <- list(log2_complete_case = cc, median_centered = median_centered, quantile_normalized = quantile_normalized)

pca_rows <- lapply(names(strategies), function(nm) {
  z <- t(strategies[[nm]])
  p <- prcomp(z, center = TRUE, scale. = TRUE)
  data.frame(
    abundance_column = rownames(p$x),
    strategy = nm,
    PC1 = p$x[, 1],
    PC2 = p$x[, 2],
    variance_PC1 = summary(p)$importance[2, 1],
    variance_PC2 = summary(p)$importance[2, 2],
    stringsAsFactors = FALSE
  )
})
pca <- do.call(rbind, pca_rows)

write.table(sample_metrics, file.path(outdir, "sample_qc_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(protein_metrics, file.path(outdir, "protein_qc_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pca, file.path(outdir, "pca_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p_missing <- ggplot(protein_metrics, aes(x = median_log2_abundance, y = missing_fraction, color = fdr_confidence)) +
  geom_point(alpha = 0.65, size = 1) +
  scale_color_manual(values = c(High = "#2166AC", Medium = "#B2182B")) +
  labs(x = "Median log2 exported abundance", y = "Missing fraction", color = "PD confidence") +
  theme_bw(base_size = 9)
ggsave(file.path(outdir, "missingness_vs_abundance.pdf"), p_missing, width = 5.2, height = 4.0)

p_pca <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(size = 1.4, alpha = 0.8) +
  facet_wrap(~ strategy, scales = "free") +
  theme_bw(base_size = 9) +
  labs(x = "PC1", y = "PC2")
ggsave(file.path(outdir, "pca_normalization_sensitivity.pdf"), p_pca, width = 8.2, height = 3.0)

manifest <- c(
  paste0("matrix_tsv\t", normalizePath(matrix_tsv)),
  paste0("proteins_all\t", nrow(x)),
  paste0("samples\t", ncol(x)),
  paste0("high_confidence_complete_case_proteins_for_pca\t", sum(complete)),
  "pca_missingness_policy\tComplete-case high-confidence proteins only; no imputation",
  "normalization_sensitivity\tlog2 exported abundance; sample median centering; quantile normalization",
  "group_labels_used\tfalse",
  paste0("R_version\t", R.version.string),
  paste0("limma_version\t", as.character(packageVersion("limma")))
)
writeLines(manifest, file.path(outdir, "analysis_manifest.tsv"))

