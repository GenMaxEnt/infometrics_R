# entropy.R
# Entropy measures for information-theoretic estimation.
#
# References:
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145.
#   Shannon, C.E. (1948). A Mathematical Theory of Communication. Bell
#     System Technical Journal, 27, 379-423.

# ---- helpers ----------------------------------------------------------------

#' Validate a probability vector
#'
#' @param p Numeric vector. Must be positive and sum to 1.
#' @param tol Numeric. Tolerance for the sum-to-one check. Default 1e-8.
#' @param name Character. Name used in error messages.
#' @return Invisibly returns \code{p} after validation.
#' @keywords internal
.check_prob <- function(p, tol = 1e-8, name = "p") {
  if (!is.numeric(p))         stop(name, " must be a numeric vector.")
  if (any(!is.finite(p)))     stop(name, " must contain only finite values.")
  if (any(p < 0))             stop(name, " must be non-negative.")
  if (abs(sum(p) - 1) > tol)  stop(name, " must sum to 1 (sum = ", round(sum(p), 8), ").")
  invisible(p)
}

# ---- Shannon entropy --------------------------------------------------------

#' Shannon Entropy
#'
#' Computes the Shannon entropy H(p) = -sum(p * log(p)), using the natural
#' logarithm (nats) by default. The convention 0 * log(0) = 0 is applied.
#'
#' @param p Numeric vector of probabilities. Must be non-negative and sum to 1.
#' @param base Positive numeric. Logarithm base. Use \code{base = 2} for bits,
#'   \code{base = exp(1)} (default) for nats, \code{base = 10} for hartleys.
#'
#' @return A single non-negative numeric value. Returns 0 for a degenerate
#'   (point-mass) distribution and \code{log(length(p))} for the uniform
#'   distribution.
#'
#' @details
#' Shannon's entropy measure satisfies three axioms: normalization (maximum
#' entropy for uniform distributions), continuity, and additivity for
#' independent sub-systems. It equals zero if and only if all probability
#' mass is concentrated on a single outcome (perfect certainty) and reaches
#' its maximum value of \code{log(K)} for the uniform distribution over K
#' outcomes. See Golan (2008, Section 3.1) for a detailed treatment.
#'
#' @references
#' Golan, A. (2008). Information and Entropy Econometrics. \emph{Foundations
#' and Trends in Econometrics}, \strong{2}(1-2), 1-145.
#'
#' Shannon, C.E. (1948). A mathematical theory of communication. \emph{Bell
#' System Technical Journal}, \strong{27}, 379-423.
#'
#' @examples
#' # Uniform distribution over 6 outcomes (maximum entropy)
#' shannon_entropy(rep(1/6, 6))
#'
#' # Non-uniform distribution (lower entropy)
#' shannon_entropy(c(0.5, 0.3, 0.2))
#'
#' # In bits (base 2)
#' shannon_entropy(rep(1/2, 2), base = 2)  # = 1 bit
#'
#' @export
shannon_entropy <- function(p, base = exp(1)) {
  .check_prob(p)
  if (!is.numeric(base) || length(base) != 1 || base <= 0)
    stop("base must be a single positive number.")
  # Apply 0 * log(0) = 0 convention
  idx <- p > 0
  h <- -sum(p[idx] * log(p[idx]))
  if (base != exp(1)) h <- h / log(base)
  h
}
