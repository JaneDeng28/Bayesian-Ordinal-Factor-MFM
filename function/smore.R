# =============================================================================
# Public R interface for SMORE
# Single-cell MethylOme Reduction and Embedding
# =============================================================================

# Load the statistical implementation from the same directory as this file.
# Resolving the path here allows users to source smore.R from any working
# directory.
.smore_source_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
  error = function(e) NA_character_
)
.smore_function_dir <- if (!is.na(.smore_source_file)) {
  dirname(.smore_source_file)
} else {
  file.path(getwd(), "function")
}
source(file.path(.smore_function_dir, "ordinal_mfm.R"), local = FALSE)
rm(.smore_source_file, .smore_function_dir)

#' Fit SMORE to ordinal methylation data
#'
#' This is the primary user-facing fitting function. It fits the Bayesian
#' ordinal factor mixture-of-finite-mixtures model implemented by
#' fit_ordinal_mfm() and returns an object with both the new SMORE class and the
#' legacy class for backward compatibility.
#'
#' @param Y Integer feature-by-cell ordinal matrix with entries in {1,...,5}.
#' @param ... Additional arguments passed to fit_ordinal_mfm().
#' @return An object of class smore_fit and ordinal_mfm_fit.
fit_smore <- function(Y, ...) {
  fit <- fit_ordinal_mfm(Y, ...)
  class(fit) <- unique(c("smore_fit", class(fit)))
  fit
}

#' Summarize posterior clustering from a SMORE fit
#'
#' @param fit Object returned by fit_smore(). Legacy ordinal_mfm_fit objects are
#'   also accepted.
#' @return Representative cell labels, posterior cluster-count summaries,
#'   posterior co-clustering probabilities, and Dahl diagnostics.
summarize_smore <- function(fit) {
  summarize_clustering(fit)
}

#' Fit SMORE from continuous methylation proportions
#'
#' The input is converted to five ordered categories before fitting. Missing
#' values must be handled before this function is called.
#'
#' @param W Numeric feature-by-cell methylation proportion matrix in [0,1].
#' @param ... Additional arguments passed to fit_smore().
#' @return An object of class smore_fit and ordinal_mfm_fit.
fit_smore_proportions <- function(W, ...) {
  fit_smore(methylation_to_ordinal(W), ...)
}
