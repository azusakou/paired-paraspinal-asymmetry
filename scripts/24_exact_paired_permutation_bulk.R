#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
input_xlsx <- if (length(args) >= 1) args[[1]] else "/Volumes/wang/paper14/ais_paraspinal_ecology/raw/bulk/GSE254300_readcount.xlsx"
edge_dir <- if (length(args) >= 2) args[[2]] else "results/tissue/gse254300_bulk"
outdir <- if (length(args) >= 3) args[[3]] else "results/tissue/gse254300_bulk/permutation_sensitivity"
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
meta$donor <- paste(meta$group, meta$donor_number, sep = "_")

y <- DGEList(counts)
filter_design <- model.matrix(~ 0 + factor(meta$group):factor(meta$side_code))
keep <- filterByExpr(y, design = filter_design, min.count = 10)
y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
logcpm <- cpm(y, log = TRUE, prior.count = 2)

donors <- unique(meta$donor)
delta <- vapply(donors, function(donor) {
  concave <- which(meta$donor == donor & meta$side_code == "A")
  convex <- which(meta$donor == donor & meta$side_code == "T")
  stopifnot(length(concave) == 1, length(convex) == 1)
  logcpm[, concave] - logcpm[, convex]
}, numeric(nrow(logcpm)))
rownames(delta) <- rownames(logcpm)
colnames(delta) <- donors
donor_group <- sub("_.*", "", donors)
stopifnot(sum(donor_group == "IS") == 5, sum(donor_group == "CS") == 5)

row_t <- function(x) {
  means <- rowMeans(x)
  sds <- apply(x, 1, sd)
  means / pmax(sds / sqrt(ncol(x)), .Machine$double.eps)
}
observed_t <- row_t(delta)
sign_grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), ncol(delta))))
extreme <- integer(nrow(delta))
for (i in seq_len(nrow(sign_grid))) {
  permuted_t <- row_t(sweep(delta, 2, sign_grid[i, ], `*`))
  extreme <- extreme + as.integer(abs(permuted_t) >= abs(observed_t) - 1e-12)
}
shared_p <- extreme / nrow(sign_grid)

observed_interaction <- rowMeans(delta[, donor_group == "IS", drop = FALSE]) -
  rowMeans(delta[, donor_group == "CS", drop = FALSE])
ais_assignments <- combn(seq_along(donors), 5)
interaction_extreme <- integer(nrow(delta))
for (i in seq_len(ncol(ais_assignments))) {
  ais_idx <- ais_assignments[, i]
  cs_idx <- setdiff(seq_along(donors), ais_idx)
  permuted <- rowMeans(delta[, ais_idx, drop = FALSE]) - rowMeans(delta[, cs_idx, drop = FALSE])
  interaction_extreme <- interaction_extreme + as.integer(abs(permuted) >= abs(observed_interaction) - 1e-12)
}
interaction_p <- interaction_extreme / ncol(ais_assignments)

result <- data.frame(
  gene_symbol = rownames(delta),
  mean_concave_minus_convex = rowMeans(delta),
  paired_t = observed_t,
  sign_flip_PValue = shared_p,
  sign_flip_FDR = p.adjust(shared_p, method = "BH"),
  AIS_minus_CS_side_effect = observed_interaction,
  group_permutation_PValue = interaction_p,
  group_permutation_FDR = p.adjust(interaction_p, method = "BH"),
  stringsAsFactors = FALSE
)
result <- result[order(result$sign_flip_PValue, -abs(result$mean_concave_minus_convex)), ]
write.table(result, file.path(outdir, "exact_permutation_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

edge_shared <- read.delim(file.path(edge_dir, "shared_concave_minus_convex.tsv"))
edge_interaction <- read.delim(file.path(edge_dir, "side_effect_interaction_AIS_minus_CS.tsv"))
comparison <- merge(
  edge_shared[, c("gene_symbol", "logFC", "FDR")],
  result[, c("gene_symbol", "mean_concave_minus_convex", "sign_flip_PValue", "sign_flip_FDR", "group_permutation_FDR")],
  by = "gene_symbol"
)
edge_interaction_fdr <- edge_interaction$FDR[match(comparison$gene_symbol, edge_interaction$gene_symbol)]
summary <- data.frame(
  metric = c(
    "tested_genes", "edgeR_shared_fdr05", "exact_sign_flip_shared_fdr05",
    "minimum_two_sided_sign_flip_p", "edgeR_shared_fdr05_exact_p_le_0.01",
    "edgeR_shared_fdr05_exact_p_le_0.05", "edgeR_and_exact_BH_fdr05_overlap",
    "edgeR_and_exact_BH_overlap_direction_concordant",
    "edgeR_interaction_fdr10", "exact_group_permutation_interaction_fdr10",
    "effect_spearman"
  ),
  value = c(
    nrow(result), sum(comparison$FDR < 0.05), sum(comparison$sign_flip_FDR < 0.05),
    min(comparison$sign_flip_PValue),
    sum(comparison$FDR < 0.05 & comparison$sign_flip_PValue <= 0.01),
    sum(comparison$FDR < 0.05 & comparison$sign_flip_PValue <= 0.05),
    sum(comparison$FDR < 0.05 & comparison$sign_flip_FDR < 0.05),
    sum(comparison$FDR < 0.05 & comparison$sign_flip_FDR < 0.05 & sign(comparison$logFC) == sign(comparison$mean_concave_minus_convex)),
    sum(edge_interaction_fdr < 0.10), sum(comparison$group_permutation_FDR < 0.10),
    cor(comparison$logFC, comparison$mean_concave_minus_convex, method = "spearman")
  )
)
write.table(summary, file.path(outdir, "permutation_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(donor = donors, group = donor_group), file.path(outdir, "donor_order.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  "shared_null\twithin-donor side labels are exchangeable; exhaustive 2^10 sign flips",
  "interaction_null\tdonor disease labels are exchangeable; exhaustive choose(10,5)=252 assignments",
  "expression\tTMM-normalized log2 CPM with prior count 2",
  "inference_boundary\texact randomization sensitivity analysis; primary count model remains robust edgeR quasi-likelihood",
  "multiplicity_boundary\ttwo-sided sign-flip p-values are discrete with minimum 2/1024; gene-level BH significance is underpowered and nominal exact p-values are used only as sensitivity evidence"
), file.path(outdir, "analysis_manifest.tsv"))

cat("Exact paired permutation sensitivity complete\n")
print(summary)
