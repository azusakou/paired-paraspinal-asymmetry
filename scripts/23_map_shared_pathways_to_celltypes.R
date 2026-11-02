#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(msigdbr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
expression_file <- if (length(args) >= 1) args[[1]] else "results/tissue/sc_reference/gene_sample_celltype_expression.tsv.gz"
camera_file <- if (length(args) >= 2) args[[2]] else "results/tissue/gse254300_bulk/camera_shared_concave_minus_convex.tsv"
outdir <- if (length(args) >= 3) args[[3]] else "results/integration/orthogonal_axes/celltype_pathways"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

expression <- read.delim(gzfile(expression_file), check.names = FALSE)
camera <- read.delim(camera_file, check.names = FALSE)
target_pathways <- camera$pathway[camera$FDR < 0.05]

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
sets <- split(
  c(hallmark$gene_symbol, reactome$gene_symbol),
  c(paste0("HALLMARK__", hallmark$gs_name), paste0("REACTOME__", reactome$gs_name))
)
sets <- lapply(sets[target_pathways], unique)

# Center each gene within a cell type across the four paired samples. This
# removes baseline differences between cell types before averaging modules.
group_key <- paste(expression$gene_symbol, expression$cell_type, sep = "__")
expression$centered_log2_cpm <- expression$log2_cpm - ave(expression$log2_cpm, group_key, FUN = mean)

scores <- list()
index <- 1L
for (pathway in names(sets)) {
  genes <- sets[[pathway]]
  subset <- expression[expression$gene_symbol %in% genes, ]
  overlap <- length(unique(subset$gene_symbol))
  if (overlap < 10) next
  split_rows <- split(subset, paste(subset$sample_id, subset$cell_type, sep = "__"))
  for (rows in split_rows) {
    scores[[index]] <- data.frame(
      pathway = pathway,
      sample_id = rows$sample_id[1], donor = rows$donor[1], side = rows$side[1],
      cell_type = rows$cell_type[1], n_genes = overlap,
      module_score = mean(rows$centered_log2_cpm),
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }
}
scores <- do.call(rbind, scores)
write.table(scores, file.path(outdir, "sample_celltype_pathway_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

paired <- merge(
  scores[scores$side == "concave", ],
  scores[scores$side == "convex", ],
  by = c("pathway", "donor", "cell_type", "n_genes"), suffixes = c("_concave", "_convex")
)
paired$concave_minus_convex <- paired$module_score_concave - paired$module_score_convex
write.table(paired, file.path(outdir, "donor_paired_celltype_pathway_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_rows <- split(paired, paste(paired$pathway, paired$cell_type, sep = "__"))
summary <- do.call(rbind, lapply(summary_rows, function(x) data.frame(
  pathway = x$pathway[1], cell_type = x$cell_type[1], n_genes = x$n_genes[1],
  mean_concave_minus_convex = mean(x$concave_minus_convex),
  donor_direction_consistent = length(unique(sign(x$concave_minus_convex))) == 1,
  donors_positive = sum(x$concave_minus_convex > 0),
  donors_negative = sum(x$concave_minus_convex < 0),
  stringsAsFactors = FALSE
)))
summary <- merge(summary, camera[, c("pathway", "Direction", "FDR")], by = "pathway", all.x = TRUE)
summary$bulk_sc_direction_concordant <- with(summary,
  (Direction == "Up" & mean_concave_minus_convex > 0) |
    (Direction == "Down" & mean_concave_minus_convex < 0)
)
write.table(summary, file.path(outdir, "celltype_pathway_side_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

top_pathways <- camera$pathway[camera$FDR < 0.01][1:min(12, sum(camera$FDR < 0.01))]
plot_dat <- summary[summary$pathway %in% top_pathways, ]
plot_dat$label <- sub("^(HALLMARK|REACTOME)__", "", plot_dat$pathway)
p <- ggplot(plot_dat, aes(cell_type, label, fill = mean_concave_minus_convex)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "Concave − convex\nmodule score")
ggsave(file.path(outdir, "shared_pathway_celltype_heatmap.pdf"), p, width = 10, height = 6.5)
ggsave(file.path(outdir, "shared_pathway_celltype_heatmap.png"), p, width = 10, height = 6.5, dpi = 220)

writeLines(c(
  "inference\tdescriptive donor-paired module effects; n=2 AIS single-cell donors; no p-values",
  "scoring\tmean gene-centered log2 CPM within each sample and cell type",
  "pathway_selection\tCAMERA FDR<0.05 in shared AIS/CS bulk side contrast",
  paste0("MSigDB_version\t", paste(unique(c(hallmark$db_version, reactome$db_version)), collapse = ";"))
), file.path(outdir, "analysis_manifest.tsv"))

cat("Mapped", length(unique(summary$pathway)), "shared pathways to", length(unique(summary$cell_type)), "cell types\n")
