# =============================================================================
# User-facing R interface for the Bayesian ordinal factor MFM model
# =============================================================================
# Typical workflow:
#   1. source("function/ordinal_mfm.R")
#   2. fit <- fit_ordinal_mfm(Y, seed = 123)
#   3. result <- summarize_clustering(fit)
#
# Matrix convention used throughout this file:
#   rows    = methylation features / genomic regions
#   columns = cells
#
# Most users need only fit_ordinal_mfm() and summarize_clustering(). The other
# functions handle compilation, input validation, or continuous-to-ordinal
# conversion. Computationally intensive sampling is implemented in mcmc.cpp.

#' Compile and load the C++ sampler
#'
#' This helper is called automatically by fit_ordinal_mfm() when the compiled
#' functions are not already available in the current R session. Users normally
#' do not need to call it directly.
#'
#' @param project_dir Top-level project directory containing function/mcmc.cpp.
#' @param rebuild Force recompilation even when an Rcpp cache is available.
#' @param quiet Suppress compiler output when TRUE.
#' @return Invisibly returns TRUE after successful compilation.

compile_ordinal_mfm <- function(project_dir = ".", rebuild = FALSE, quiet = TRUE) {
  # Resolve the C++ source relative to the repository, not the user's home path.
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  cpp_file <- file.path(project_dir, "function", "mcmc.cpp")
  if (!file.exists(cpp_file)) stop("Cannot find function/mcmc.cpp under project_dir.", call. = FALSE)
  if (!requireNamespace("Rcpp", quietly = TRUE) ||
      !requireNamespace("RcppArmadillo", quietly = TRUE)) {
    stop("Install Rcpp and RcppArmadillo before fitting the model.", call. = FALSE)
  }

  # Compilation artifacts are kept in R's temporary directory so they are not
  # added to the repository or confused with scientific output.
  cache_dir <- file.path(tempdir(), "ordinal_mfm_cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  old_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  temporary_makevars <- NULL
  on.exit({
    if (is.na(old_makevars)) Sys.unsetenv("R_MAKEVARS_USER") else Sys.setenv(R_MAKEVARS_USER = old_makevars)
    if (!is.null(temporary_makevars)) unlink(temporary_makevars)
  }, add = TRUE)
  # Compatibility fallback for macOS R installations whose compiler settings
  # still point to /opt/gfortran although the runtime is bundled in R.framework.
  # The user's original R_MAKEVARS_USER value is restored on function exit.
  if (identical(Sys.info()[["sysname"]], "Darwin") &&
      !dir.exists("/opt/gfortran") &&
      file.exists(file.path(R.home("lib"), "libgfortran.5.dylib"))) {
    temporary_makevars <- tempfile("ordinal-mfm-Makevars-")
    writeLines(sprintf("FLIBS = %s %s",
                       file.path(R.home("lib"), "libgfortran.5.dylib"),
                       file.path(R.home("lib"), "libquadmath.0.dylib")),
               temporary_makevars)
    Sys.setenv(R_MAKEVARS_USER = temporary_makevars)
  }
  Rcpp::sourceCpp(cpp_file, cacheDir = cache_dir, rebuild = rebuild,
                  showOutput = !quiet, verbose = FALSE)
  invisible(TRUE)
}

#' Convert methylation proportions to five ordered categories
#'
#' The cut intervals are [0,0.2], (0.2,0.4], (0.4,0.6], (0.6,0.8], and
#' (0.8,1]. This matches the ordinal representation described in the paper.
#'
#' @param W Numeric feature-by-cell matrix of methylation proportions in [0,1].
#' @return Integer matrix with the same dimensions/dimnames and values 1,...,5.
methylation_to_ordinal <- function(W) {
  W <- as.matrix(W)
  if (!is.numeric(W) || anyNA(W) || any(!is.finite(W))) {
    stop("W must be a numeric matrix without NA, NaN, or Inf.", call. = FALSE)
  }
  if (any(W < 0 | W > 1)) stop("W must contain methylation proportions in [0, 1].", call. = FALSE)
  matrix(findInterval(W, c(0.2, 0.4, 0.6, 0.8), left.open = TRUE) + 1L,
         nrow = nrow(W), ncol = ncol(W), dimnames = dimnames(W))
}

# Internal input validator. It is intentionally kept separate from the sampler
# so invalid data fail before C++ compilation or a long MCMC run begins.
validate_ordinal_matrix <- function(Y) {
  Y <- as.matrix(Y)
  if (anyNA(Y) || any(!is.finite(Y)) || !is.numeric(Y)) {
    stop("Y must be a numeric/integer matrix without missing or non-finite values.", call. = FALSE)
  }
  if (nrow(Y) < 1L || ncol(Y) < 2L) stop("Y must have features in rows and at least two cells in columns.", call. = FALSE)
  if (any(Y != as.integer(Y)) || any(!Y %in% 1:5)) {
    stop("Every entry of Y must be one of the ordinal categories 1, 2, 3, 4, or 5.", call. = FALSE)
  }
  storage.mode(Y) <- "integer"
  Y
}

#' Fit the Bayesian ordinal factor MFM model
#'
#' Four independent chains are run sequentially by default. Each chain retains
#' n_iter - burn_in posterior draws. The number of clusters is inferred by the
#' MFM prior; q specifies the latent factor dimension and is not a cluster count.
#'
#' @param Y Integer feature-by-cell ordinal matrix with entries in {1,...,5}.
#' @param q Latent factor dimension (paper default: 4).
#' @param chains Number of independent MCMC chains (default: 4).
#' @param n_iter Total iterations per chain (default: 2000).
#' @param burn_in Iterations discarded at the start of each chain (default: 1000).
#' @param initial_clusters Cluster count used only to initialize the sampler.
#' @param sigma_B Prior standard deviation for free factor-loading elements.
#' @param alpha_mfm MFM Dirichlet concentration parameter.
#' @param max_clusters Computational upper bound on occupied clusters.
#' @param kappa0 Normal-inverse-Wishart mean-precision hyperparameter.
#' @param nu0_add Sets inverse-Wishart degrees of freedom to q + nu0_add.
#' @param lambda_pois Poisson prior mean parameter for the component count minus one.
#' @param S0_scale_mult Multiplier applied to the data-informed S0 initialization.
#' @param K_sum_max Truncation limit used when evaluating the MFM V_n terms.
#' @param loglik_every Compute saved log likelihood at this interval.
#' @param seed Optional base seed. Chain c uses seed + c - 1.
#' @param verbose Print chain and sampler progress.
#' @param print_every Sampler progress-printing interval.
#' @param project_dir Top-level directory containing function/mcmc.cpp.
#' @return An object of class ordinal_mfm_fit containing raw results for every
#'   chain, input dimensions, and the settings/seeds used for the analysis.
fit_ordinal_mfm <- function(
    Y, q = 4L, chains = 4L, n_iter = 2000L, burn_in = 1000L, initial_clusters = 10L,
    sigma_B = 5, alpha_mfm = 1, max_clusters = 40L, kappa0 = 0.1,
    nu0_add = 2, lambda_pois = 9, S0_scale_mult = 1,
    K_sum_max = 200L, loglik_every = 1L, seed = NULL,
    verbose = TRUE, print_every = 100L, project_dir = ".") {
  Y <- validate_ordinal_matrix(Y)
  if (q < 1L || q > nrow(Y) || q != as.integer(q)) {
    stop("q must be an integer between 1 and nrow(Y).", call. = FALSE)
  }
  if (chains < 1L || chains != as.integer(chains)) {
    stop("chains must be a positive integer.", call. = FALSE)
  }
  # Compile lazily: sourcing this R file does not trigger a compiler or MCMC run.
  if (!exists("ordinal_factor_mfm_mcmc", mode = "function")) {
    compile_ordinal_mfm(project_dir = project_dir)
  }
  # Use reproducible consecutive seeds when the user supplies a base seed.
  if (is.null(seed)) {
    chain_seeds <- sample.int(.Machine$integer.max, chains)
  } else {
    if (length(seed) != 1L || !is.finite(seed)) stop("seed must be one finite integer.", call. = FALSE)
    chain_seeds <- as.integer(seed) + seq_len(chains) - 1L
  }
  # Chains are run sequentially to avoid silently multiplying peak memory use
  # for large feature-by-cell matrices.
  raw_fits <- vector("list", chains)
  for (chain_id in seq_len(chains)) {
    if (verbose) message(sprintf("Starting MCMC chain %d of %d (seed %d)",
                                 chain_id, chains, chain_seeds[chain_id]))
    set.seed(chain_seeds[chain_id])
    raw_fits[[chain_id]] <- ordinal_factor_mfm_mcmc(
      Y = Y, K = as.integer(q), n_iter = as.integer(n_iter),
      n_burn = as.integer(burn_in), K_init = as.integer(initial_clusters),
      sigma_B = sigma_B, alpha_mfm = alpha_mfm,
      K_max_cluster = as.integer(max_clusters), kappa0 = kappa0,
      nu0_add = nu0_add, lambda_pois = lambda_pois,
      S0_scale_mult = S0_scale_mult, K_sum_max = as.integer(K_sum_max),
      loglik_every = as.integer(loglik_every), verbose = verbose,
      print_every = as.integer(print_every)
    )
  }
  names(raw_fits) <- paste0("chain_", seq_len(chains))
  structure(list(
    chains = raw_fits,
    input = list(n_features = nrow(Y), n_cells = ncol(Y), q = as.integer(q)),
    settings = list(chains = as.integer(chains), n_iter = as.integer(n_iter),
                    burn_in = as.integer(burn_in), chain_seeds = chain_seeds,
                    initial_clusters = as.integer(initial_clusters), sigma_B = sigma_B,
                    alpha_mfm = alpha_mfm, max_clusters = as.integer(max_clusters),
                    kappa0 = kappa0, nu0_add = nu0_add, lambda_pois = lambda_pois,
                    S0_scale_mult = S0_scale_mult, K_sum_max = as.integer(K_sum_max))
  ), class = "ordinal_mfm_fit")
}

#' Summarize posterior clustering across all chains
#'
#' Post-burn-in allocation draws from all chains are combined. The function
#' computes the posterior similarity matrix and selects Dahl's least-squares
#' representative partition. Per-chain cluster-count traces are retained for
#' convergence and sensitivity checks.
#'
#' @param fit Object returned by fit_ordinal_mfm().
#' @return A list containing cell labels, inferred cluster counts, posterior
#'   co-clustering probabilities, per-chain traces/modes, and Dahl diagnostics.
summarize_clustering <- function(fit) {
  if (!inherits(fit, "ordinal_mfm_fit")) stop("fit must be returned by fit_ordinal_mfm().", call. = FALSE)
  if (!length(fit$chains)) stop("The fit contains no MCMC chains.", call. = FALSE)
  if (!exists("dahl_partition_cpp", mode = "function")) {
    stop("Compiled Dahl function is unavailable; rerun fit_ordinal_mfm().", call. = FALSE)
  }
  # Each z_trace has cells in rows and retained iterations in columns.
  draws_by_chain <- lapply(fit$chains, function(raw) {
    if (raw$n_save < 1L) stop("A chain contains no post-burn-in samples.", call. = FALSE)
    raw$z_trace[, seq_len(raw$n_save), drop = FALSE]
  })
  # Combining allocation draws is valid for the co-clustering representation:
  # it is invariant to arbitrary numeric cluster labels across chains.
  draws <- do.call(cbind, draws_by_chain)
  dahl <- dahl_partition_cpp(draws)
  labels <- as.integer(dahl$z_dahl)
  ncluster_by_chain <- lapply(fit$chains, function(raw) raw$nclust_trace[seq_len(raw$n_save)])
  ncluster_draws <- unlist(ncluster_by_chain, use.names = FALSE)
  list(
    cluster = labels,
    n_clusters = length(unique(labels)),
    posterior_mode_n_clusters = as.integer(names(which.max(table(ncluster_draws)))),
    posterior_similarity = dahl$psm,
    ncluster_trace = ncluster_draws,
    ncluster_trace_by_chain = ncluster_by_chain,
    posterior_mode_by_chain = vapply(ncluster_by_chain, function(x) {
      as.integer(names(which.max(table(x))))
    }, integer(1)),
    dahl_loss = dahl$loss,
    selected_draw = dahl$best_m
  )
}

#' Fit the model directly from methylation proportions
#'
#' Convenience wrapper equivalent to calling methylation_to_ordinal(W) and then
#' fit_ordinal_mfm(). Missing-value handling must be completed beforehand.
#'
#' @param W Numeric feature-by-cell methylation proportion matrix in [0,1].
#' @param ... Additional arguments passed to fit_ordinal_mfm().
#' @return An ordinal_mfm_fit object.
fit_methylation_mfm <- function(W, ...) {
  fit_ordinal_mfm(methylation_to_ordinal(W), ...)
}

# Backward-compatible alias from an earlier draft. New analyses should call
# fit_ordinal_mfm() for ordinal Y or fit_methylation_mfm() for continuous W.
run_analysis <- function(W, K = 10L, n_burn = 1000L, ...) {
  fit_methylation_mfm(W, q = K, burn_in = n_burn, ...)
}
