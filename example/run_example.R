#!/usr/bin/env Rscript

# =============================================================================
# End-to-end example: fit the paper's model to simulated balanced/strong data
# =============================================================================
# Run from a terminal with:
#   Rscript example/run_example.R
#
# This is a full-size example rather than a quick unit test. The default call
# runs four sequential chains of 2,000 iterations each on a 5,000 x 1,000
# matrix. Runtime and memory requirements will therefore be substantial.

# Balanced, strong-signal model-based simulation from the paper:
# p = 5000 features, n = 1000 cells, 5 balanced true clusters.
# Resolve paths from this script's location so the example works regardless of
# the terminal's current working directory.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
project_dir <- normalizePath(file.path(dirname(normalizePath(script_file)), ".."))

# Load the public R functions. The C++ backend compiles automatically at the
# first fit_smore() call in a new R session.
source(file.path(project_dir, "function", "smore.R"))

# The RDS object contains:
#   Y          5,000 x 1,000 ordinal input matrix (values 1,...,5)
#   W          corresponding continuous methylation proportions
#   true_label simulated labels, used only to evaluate this example
#   metadata   complete data-generation settings
example_data <- readRDS(file.path(
  project_dir, "data", "example_data_n1000_p5000.rds"
))
# The proposed model uses Y; true_label is never supplied to the sampler.
Y <- example_data$Y

# Paper defaults are used automatically:
# q = 4, chains = 4, n_iter = 2000, burn_in = 1000,
# initial_clusters = 10, sigma_B = 5, alpha_mfm = 1,
# max_clusters = 40, kappa0 = 0.1, nu0 = q + 2, lambda_pois = 9.
# Consecutive chain seeds are 20260820, 20260821, 20260822, and 20260823.
fit <- fit_smore(Y, seed = 20260820, project_dir = project_dir)

# Combine retained allocation draws from all chains, build the posterior
# similarity matrix, and select Dahl's least-squares representative partition.
clustering <- summarize_smore(fit)

# ARI is available here because simulated truth is included with the example.
ari <- mclust::adjustedRandIndex(example_data$true_label, clustering$cluster)
metrics <- data.frame(
  metric = c("dahl_clusters", "posterior_mode_clusters", "adjusted_rand_index"),
  value = c(clustering$n_clusters, clustering$posterior_mode_n_clusters, ari)
)
assignments <- data.frame(
  cell = colnames(Y),
  true_cluster = example_data$true_label,
  estimated_cluster = clustering$cluster
)

# Save human-readable summaries as CSV and the complete fit as RDS. The RDS file
# contains all four raw chains plus the combined posterior clustering summary.
write.csv(metrics, file.path(project_dir, "example", "example_metrics.csv"), row.names = FALSE)
write.csv(assignments, file.path(project_dir, "example", "example_cluster_assignments.csv"),
          row.names = FALSE)
saveRDS(list(fit = fit, clustering = clustering),
        file.path(project_dir, "example", "example_results.rds"), compress = "xz")

# Print overall metrics followed by the posterior modal cluster count from each
# chain. Large disagreements across chains should be investigated before using
# the representative partition for scientific interpretation.
print(metrics, row.names = FALSE)
print(clustering$posterior_mode_by_chain)
message("Results saved in: ", file.path(project_dir, "example"))
