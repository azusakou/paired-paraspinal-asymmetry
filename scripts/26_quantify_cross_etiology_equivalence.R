#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
bulk_dir <- if (length(args) >= 1) args[[1]] else "results/tissue/gse254300_bulk"
outdir <- if (length(args) >= 2) args[[2]] else file.path(bulk_dir, "cross_etiology_equivalence")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(stem, prefix) {
  x <- read.delim(file.path(bulk_dir, paste0(stem, ".tsv")), check.names = FALSE)
  needed <- c(
    "gene_symbol", "logFC", "standard_error", "df_total",
    "ci90_low", "ci90_high", "logCPM", "PValue", "FDR"
  )
  stopifnot(all(needed %in% names(x)))
  x <- x[, needed]
  names(x)[-1] <- paste0(prefix, "_", names(x)[-1])
  x
}

shared <- read_result("shared_concave_minus_convex", "shared")
ais <- read_result("AIS_concave_minus_convex", "ais")
cs <- read_result("CS_concave_minus_convex", "cs")
interaction <- read_result("side_effect_interaction_AIS_minus_CS", "interaction")
result <- Reduce(function(x, y) merge(x, y, by = "gene_symbol"), list(shared, ais, cs, interaction))

result$same_direction <- sign(result$ais_logFC) == sign(result$cs_logFC)
result$shared_fdr05 <- result$shared_FDR < 0.05

# Inverse-variance fixed-effect summaries and a descriptive two-stratum I2.
# The formal disease-by-side QL interaction remains the primary heterogeneity test.
wa <- 1 / result$ais_standard_error^2
wc <- 1 / result$cs_standard_error^2
result$fixed_effect_logFC <- (wa * result$ais_logFC + wc * result$cs_logFC) / (wa + wc)
result$fixed_effect_standard_error <- sqrt(1 / (wa + wc))
result$heterogeneity_Q <- wa * (result$ais_logFC - result$fixed_effect_logFC)^2 +
  wc * (result$cs_logFC - result$fixed_effect_logFC)^2
result$heterogeneity_I2 <- ifelse(
  result$heterogeneity_Q > 0,
  pmax(0, (result$heterogeneity_Q - 1) / result$heterogeneity_Q), 0
)

tost <- function(estimate, standard_error, df, margin) {
  lower_p <- pt((estimate + margin) / standard_error, df = df, lower.tail = FALSE)
  upper_p <- pt((estimate - margin) / standard_error, df = df, lower.tail = TRUE)
  pmax(lower_p, upper_p)
}

margins <- c(0.25, 0.50, 1.00)
summary_rows <- list()
for (i in seq_along(margins)) {
  margin <- margins[i]
  suffix <- gsub("\\.", "p", format(margin, nsmall = 2))
  p_col <- paste0("equivalence_p_margin_", suffix)
  fdr_col <- paste0("equivalence_fdr_margin_", suffix)
  result[[p_col]] <- tost(
    result$interaction_logFC, result$interaction_standard_error,
    result$interaction_df_total, margin
  )
  result[[fdr_col]] <- p.adjust(result[[p_col]], method = "BH")
  equivalent <- result[[fdr_col]] < 0.05
  summary_rows[[i]] <- data.frame(
    equivalence_margin_log2 = margin,
    maximum_fold_difference = 2^margin,
    all_genes_equivalent_fdr05 = sum(equivalent, na.rm = TRUE),
    shared_hits = sum(result$shared_fdr05),
    shared_hits_equivalent_fdr05 = sum(result$shared_fdr05 & equivalent, na.rm = TRUE),
    shared_hits_equivalent_fraction = mean(equivalent[result$shared_fdr05], na.rm = TRUE),
    shared_hits_ci90_inside_margin = sum(
      result$shared_fdr05 & result$interaction_ci90_low > -margin &
        result$interaction_ci90_high < margin,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

result$strictly_conserved_margin_0p50 <- with(
  result,
  shared_fdr05 & same_direction & equivalence_fdr_margin_0p50 < 0.05
)
result$conservation_class <- "not_shared_fdr05"
result$conservation_class[result$shared_fdr05] <- "shared_but_equivalence_inconclusive"
result$conservation_class[result$strictly_conserved_margin_0p50] <- "strictly_conserved_margin_0p50"

result <- result[order(result$strictly_conserved_margin_0p50, result$shared_FDR), ]
write.table(result, gzfile(file.path(outdir, "gene_cross_etiology_equivalence.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
equivalence_summary <- do.call(rbind, summary_rows)
write.table(equivalence_summary, file.path(outdir, "equivalence_margin_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

overall <- data.frame(
  metric = c(
    "tested_genes", "shared_fdr05", "shared_same_direction",
    "interaction_fdr10", "strictly_conserved_margin_0p50",
    "strictly_conserved_fraction_of_shared", "AIS_CS_effect_spearman",
    "shared_median_I2", "shared_fraction_I2_below_0p25"
  ),
  value = c(
    nrow(result), sum(result$shared_fdr05),
    sum(result$shared_fdr05 & result$same_direction),
    sum(result$interaction_FDR < 0.10),
    sum(result$strictly_conserved_margin_0p50),
    mean(result$strictly_conserved_margin_0p50[result$shared_fdr05]),
    cor(result$ais_logFC, result$cs_logFC, method = "spearman"),
    median(result$heterogeneity_I2[result$shared_fdr05], na.rm = TRUE),
    mean(result$heterogeneity_I2[result$shared_fdr05] < 0.25, na.rm = TRUE)
  )
)
write.table(overall, file.path(outdir, "equivalence_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_dat <- equivalence_summary
plot_dat$label <- sprintf("+/- %.2f log2 (<= %.2f-fold)", plot_dat$equivalence_margin_log2, plot_dat$maximum_fold_difference)
p_margin <- ggplot(plot_dat, aes(label, shared_hits_equivalent_fraction)) +
  geom_col(width = 0.65, fill = "#2166AC") +
  geom_text(aes(label = sprintf("%d / %d", shared_hits_equivalent_fdr05, shared_hits)), vjust = -0.4, size = 3.2) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05), expand = expansion(mult = c(0, 0))) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(x = "Equivalence margin (retrospective sensitivity analysis)", y = "Shared hits with FDR-controlled equivalence")
ggsave(file.path(outdir, "equivalence_margin_sensitivity.pdf"), p_margin, width = 5.5, height = 4.2)
ggsave(file.path(outdir, "equivalence_margin_sensitivity.png"), p_margin, width = 5.5, height = 4.2, dpi = 240)

p_scatter <- ggplot(result, aes(ais_logFC, cs_logFC, color = conservation_class)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.35) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c(
    not_shared_fdr05 = "grey80",
    shared_but_equivalence_inconclusive = "#92C5DE",
    strictly_conserved_margin_0p50 = "#B2182B"
  ), labels = c(
    not_shared_fdr05 = "Not shared",
    shared_but_equivalence_inconclusive = "Shared; equivalence inconclusive",
    strictly_conserved_margin_0p50 = "Equivalent (+/-0.50 log2)"
  )) +
  coord_equal() + theme_bw(base_size = 9) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(
    x = "AIS concave - convex log2 fold change",
    y = "Congenital scoliosis concave - convex log2 fold change",
    color = NULL
  )
ggsave(file.path(outdir, "strict_conservation_scatter.pdf"), p_scatter, width = 7.2, height = 5.8)
ggsave(file.path(outdir, "strict_conservation_scatter.png"), p_scatter, width = 7.2, height = 5.8, dpi = 240)

writeLines(c(
  "primary_question\twhether the AIS-minus-CS side interaction is practically equivalent within a prespecified margin",
  "primary_margin\tplus/minus 0.50 log2 fold change (maximum 1.414-fold difference between etiologies); declared for this retrospective equivalence analysis but not preregistered",
  "sensitivity_margins\tplus/minus 0.25 and plus/minus 1.00 log2 fold change",
  "equivalence_test\ttwo one-sided tests using edgeR quasi-likelihood contrast estimate, standard error and total degrees of freedom",
  "multiplicity\tBenjamini-Hochberg across the full retained gene universe separately at each margin",
  "strict_conservation\tshared-side FDR<0.05; same AIS/CS direction; interaction equivalence FDR<0.05 at margin 0.50",
  "boundary\tequivalence margin is an analysis-scale smallest effect size of interest, not a clinically validated threshold",
  "heterogeneity\tformal edgeR interaction is primary; inverse-variance Q and I2 are descriptive two-stratum summaries"
), file.path(outdir, "analysis_manifest.tsv"))

cat("Cross-etiology equivalence analysis complete\n")
print(overall)
print(equivalence_summary)
