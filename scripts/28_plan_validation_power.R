#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[[1]] else "results/design/prospective_validation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Paired side-effect validation. Cohen dz is the mean paired difference divided
# by the SD of paired differences. Alpha 0.025 allows two co-primary functional
# endpoints; final multiplicity must be frozen in the prospective SAP.
paired_grid <- expand.grid(
  standardized_paired_effect = c(0.40, 0.50, 0.60, 0.80, 1.00),
  target_power = c(0.80, 0.90),
  alpha_two_sided = 0.025,
  stringsAsFactors = FALSE
)
paired_grid$complete_pairs_required <- mapply(function(effect, power, alpha) {
  ceiling(power.t.test(
    delta = effect, sd = 1, sig.level = alpha, power = power,
    type = "paired", alternative = "two.sided"
  )$n)
}, paired_grid$standardized_paired_effect, paired_grid$target_power, paired_grid$alpha_two_sided)
paired_grid$recruit_with_15pct_attrition <- ceiling(paired_grid$complete_pairs_required / 0.85)
write.table(paired_grid, file.path(outdir, "paired_effect_power_scenarios.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Two-etiology equivalence of paired side effects. Simulation assumes equal
# within-etiology SD after converting the endpoint to standardized paired
# differences, true interaction zero, and two independent etiology groups.
simulate_equivalence_power <- function(n, margin, simulations = 50000L, alpha = 0.05) {
  df0 <- n - 1
  var1 <- rchisq(simulations, df0) / df0
  var2 <- rchisq(simulations, df0) / df0
  difference <- rnorm(simulations, mean = 0, sd = sqrt(2 / n))
  standard_error <- sqrt(var1 / n + var2 / n)
  welch_df <- (var1 / n + var2 / n)^2 /
    ((var1 / n)^2 / df0 + (var2 / n)^2 / df0)
  p_lower <- pt((difference + margin) / standard_error, df = welch_df, lower.tail = FALSE)
  p_upper <- pt((difference - margin) / standard_error, df = welch_df, lower.tail = TRUE)
  mean(p_lower < alpha & p_upper < alpha)
}

set.seed(14016)
equivalence_grid <- expand.grid(
  n_per_etiology = seq(15, 100),
  standardized_equivalence_margin = c(0.50, 0.75),
  stringsAsFactors = FALSE
)
equivalence_grid$power_if_true_interaction_zero <- mapply(
  simulate_equivalence_power,
  equivalence_grid$n_per_etiology,
  equivalence_grid$standardized_equivalence_margin
)
write.table(equivalence_grid, file.path(outdir, "etiology_equivalence_power_simulation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

targets <- do.call(rbind, lapply(split(equivalence_grid, equivalence_grid$standardized_equivalence_margin), function(x) {
  do.call(rbind, lapply(c(0.80, 0.90), function(target_power) {
    eligible <- x[x$power_if_true_interaction_zero >= target_power, ]
    if (!nrow(eligible)) return(data.frame(
      standardized_equivalence_margin = x$standardized_equivalence_margin[1],
      target_power = target_power, complete_donors_per_etiology = NA,
      recruit_per_etiology_with_15pct_attrition = NA
    ))
    n <- min(eligible$n_per_etiology)
    data.frame(
      standardized_equivalence_margin = x$standardized_equivalence_margin[1],
      target_power = target_power,
      complete_donors_per_etiology = n,
      recruit_per_etiology_with_15pct_attrition = ceiling(n / 0.85)
    )
  }))
}))
write.table(targets, file.path(outdir, "etiology_equivalence_sample_targets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  "paired_primary_alpha\t0.025 two-sided for scenario planning with two co-primary endpoints",
  "paired_effect_scale\tCohen dz based on donor-level concave-minus-convex differences",
  "equivalence_alpha\t0.05 for both one-sided tests",
  "equivalence_assumption\ttrue standardized etiology interaction zero; equal SD; independent etiology groups",
  "equivalence_simulations\t50000 per sample-size and margin scenario",
  "attrition_inflation\t15 percent",
  "boundary\tscenario analysis only; replace standardized assumptions with blinded pilot variance and a scientifically justified raw-scale margin before preregistration"
), file.path(outdir, "analysis_manifest.tsv"))

cat("Prospective validation power scenarios complete\n")
print(paired_grid)
print(targets)
