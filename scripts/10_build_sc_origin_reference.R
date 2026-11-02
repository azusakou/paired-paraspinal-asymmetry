#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1) args[[1]] else "/Volumes/wang/paper15/ais_mesenchymal_axis/derived/prjna722100/qc_annotation/prjna722100_qc_annotated_sce.rds"
outdir <- if (length(args) >= 2) args[[2]] else "results/tissue/sc_reference"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(input_rds)
required_cd <- c("sample_id", "donor", "side", "coarse_label", "scDblFinder.class")
stopifnot(all(required_cd %in% colnames(colData(sce))))
stopifnot("counts" %in% assayNames(sce))

cd <- as.data.frame(colData(sce))
rd <- as.data.frame(rowData(sce))
counts <- assay(sce, "counts")
stopifnot(inherits(counts, "sparseMatrix"))
stopifnot(all(cd$scDblFinder.class == "singlet"))

clean_side <- tolower(as.character(cd$side))
clean_side[clean_side == "concave"] <- "concave"
clean_side[clean_side == "convex"] <- "convex"
stopifnot(all(clean_side %in% c("concave", "convex")))
cd$side <- clean_side

aggregate_sparse <- function(mat, group) {
  group <- factor(group)
  design <- sparse.model.matrix(~ 0 + group)
  colnames(design) <- sub("^group", "", colnames(design))
  list(matrix = mat %*% design, n = table(group))
}

feature <- data.frame(
  ensembl_id = if ("ensembl_id" %in% colnames(rd)) rd$ensembl_id else rownames(sce),
  gene_symbol = if ("symbol" %in% colnames(rd)) rd$symbol else rownames(sce),
  stringsAsFactors = FALSE
)

celltype <- factor(cd$coarse_label)
agg_ct <- aggregate_sparse(counts, celltype)
det_ct <- aggregate_sparse(counts > 0, celltype)
ct_cpm <- edgeR::cpm(agg_ct$matrix, log = FALSE, prior.count = 0)
ct_prev <- sweep(as.matrix(det_ct$matrix), 2, as.numeric(agg_ct$n), "/")

max_cpm <- apply(ct_cpm, 1, max)
sum_cpm <- rowSums(ct_cpm)
dominant_idx <- max.col(ct_cpm, ties.method = "first")
specificity_fraction <- ifelse(sum_cpm > 0, max_cpm / sum_cpm, NA_real_)
scaled <- sweep(ct_cpm, 1, pmax(max_cpm, .Machine$double.eps), "/")
tau <- rowSums(1 - scaled) / (ncol(ct_cpm) - 1)
tau[max_cpm == 0] <- NA_real_

gene_specificity <- cbind(
  feature,
  dominant_cell_type = colnames(ct_cpm)[dominant_idx],
  dominant_cpm = max_cpm,
  specificity_fraction = specificity_fraction,
  tau = tau,
  pooled_prevalence_in_dominant_type = ct_prev[cbind(seq_len(nrow(ct_prev)), dominant_idx)]
)

for (j in seq_len(ncol(ct_cpm))) {
  nm <- make.names(colnames(ct_cpm)[j])
  gene_specificity[[paste0("cpm_", nm)]] <- ct_cpm[, j]
  gene_specificity[[paste0("prevalence_", nm)]] <- ct_prev[, j]
}

donor_type <- interaction(cd$donor, cd$coarse_label, sep = "__", drop = TRUE)
agg_dt <- aggregate_sparse(counts, donor_type)
dt_cpm <- edgeR::cpm(agg_dt$matrix, log = FALSE, prior.count = 0)
donors <- sort(unique(as.character(cd$donor)))
dominant_by_donor <- lapply(donors, function(d) {
  cols <- grep(paste0("^", d, "__"), colnames(dt_cpm))
  types <- sub(paste0("^", d, "__"), "", colnames(dt_cpm)[cols])
  types[max.col(dt_cpm[, cols, drop = FALSE], ties.method = "first")]
})
names(dominant_by_donor) <- donors
for (d in donors) gene_specificity[[paste0("dominant_type_", d)]] <- dominant_by_donor[[d]]
gene_specificity$dominant_type_consistent_across_donors <- Reduce(`==`, dominant_by_donor)

sample_type <- interaction(cd$sample_id, cd$coarse_label, sep = "__", drop = TRUE)
agg_st <- aggregate_sparse(counts, sample_type)
det_st <- aggregate_sparse(counts > 0, sample_type)
st_cpm <- edgeR::cpm(agg_st$matrix, log = TRUE, prior.count = 0.5)
st_prev <- sweep(as.matrix(det_st$matrix), 2, as.numeric(agg_st$n), "/")

sample_meta <- do.call(rbind, strsplit(colnames(st_cpm), "__", fixed = TRUE))
sample_meta <- data.frame(sample_id = sample_meta[, 1], coarse_label = sample_meta[, 2], stringsAsFactors = FALSE)
sample_lookup <- unique(cd[, c("sample_id", "donor", "side")])
sample_meta <- merge(sample_meta, sample_lookup, by = "sample_id", sort = FALSE)
sample_meta$matrix_col <- match(paste(sample_meta$sample_id, sample_meta$coarse_label, sep = "__"), colnames(st_cpm))
sample_meta <- sample_meta[order(sample_meta$coarse_label, sample_meta$donor, sample_meta$side), ]

side_rows <- vector("list", length(levels(celltype)))
names(side_rows) <- levels(celltype)
for (ct in levels(celltype)) {
  sm <- sample_meta[sample_meta$coarse_label == ct, ]
  stopifnot(all(table(sm$donor) == 2))
  deltas <- sapply(donors, function(d) {
    sub <- sm[sm$donor == d, ]
    concave_col <- sub$matrix_col[sub$side == "concave"]
    convex_col <- sub$matrix_col[sub$side == "convex"]
    st_cpm[, concave_col] - st_cpm[, convex_col]
  })
  if (is.null(dim(deltas))) deltas <- matrix(deltas, ncol = 1)
  colnames(deltas) <- donors
  out <- data.frame(
    feature,
    cell_type = ct,
    mean_log2cpm_concave_minus_convex = rowMeans(deltas),
    donors_positive = rowSums(deltas > 0),
    donors_negative = rowSums(deltas < 0),
    direction_consistent = apply(deltas, 1, function(z) all(z > 0) || all(z < 0)),
    stringsAsFactors = FALSE
  )
  for (d in donors) out[[paste0("delta_", d)]] <- deltas[, d]
  side_rows[[ct]] <- out
}
side_effects <- do.call(rbind, side_rows)
rownames(side_effects) <- NULL

prevalence_long <- do.call(rbind, lapply(seq_len(ncol(st_prev)), function(j) {
  data.frame(
    feature,
    sample_id = sample_meta$sample_id[match(j, sample_meta$matrix_col)],
    donor = sample_meta$donor[match(j, sample_meta$matrix_col)],
    side = sample_meta$side[match(j, sample_meta$matrix_col)],
    cell_type = sample_meta$coarse_label[match(j, sample_meta$matrix_col)],
    n_cells = as.integer(agg_st$n[j]),
    prevalence = st_prev[, j],
    log2_cpm = st_cpm[, j],
    stringsAsFactors = FALSE
  )
}))

qc_summary <- data.frame(
  sample_id = names(table(cd$sample_id)),
  donor = vapply(names(table(cd$sample_id)), function(s) unique(cd$donor[cd$sample_id == s]), character(1)),
  side = vapply(names(table(cd$sample_id)), function(s) unique(cd$side[cd$sample_id == s]), character(1)),
  n_cells = as.integer(table(cd$sample_id)),
  median_umi = vapply(names(table(cd$sample_id)), function(s) median(cd$sum[cd$sample_id == s]), numeric(1)),
  median_genes = vapply(names(table(cd$sample_id)), function(s) median(cd$detected[cd$sample_id == s]), numeric(1)),
  median_mito_percent = vapply(names(table(cd$sample_id)), function(s) median(cd$subsets_Mito_percent[cd$sample_id == s]), numeric(1)),
  stringsAsFactors = FALSE
)

composition <- as.data.frame(table(sample_id = cd$sample_id, donor = cd$donor, side = cd$side, cell_type = cd$coarse_label))
composition <- composition[composition$Freq > 0, ]
composition$fraction <- ave(composition$Freq, composition$sample_id, FUN = function(z) z / sum(z))

write.table(gene_specificity, file.path(outdir, "gene_celltype_specificity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(side_effects, file.path(outdir, "gene_celltype_side_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(prevalence_long, file.path(outdir, "gene_sample_celltype_expression.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qc_summary, file.path(outdir, "qc_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(composition, file.path(outdir, "celltype_composition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot(composition, aes(x = sample_id, y = fraction, fill = cell_type)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Cell fraction", fill = "Coarse cell type") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(outdir, "celltype_composition.pdf"), p, width = 6.2, height = 3.6)
ggsave(file.path(outdir, "celltype_composition.png"), p, width = 6.2, height = 3.6, dpi = 300)

manifest <- c(
  paste0("input_rds\t", normalizePath(input_rds)),
  paste0("genes\t", nrow(sce)),
  paste0("cells\t", ncol(sce)),
  paste0("donors\t", length(donors)),
  paste0("samples\t", length(unique(cd$sample_id))),
  "statistical_boundary\tSide effects are descriptive donor-paired log2-CPM differences; no cell-level p-values are reported.",
  paste0("R_version\t", R.version.string),
  paste0("SingleCellExperiment_version\t", as.character(packageVersion("SingleCellExperiment"))),
  paste0("edgeR_version\t", as.character(packageVersion("edgeR")))
)
writeLines(manifest, file.path(outdir, "analysis_manifest.tsv"))
