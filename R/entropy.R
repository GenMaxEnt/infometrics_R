# entropy.R
# Entropy and divergence measures for information-theoretic estimation.
#
# References:
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145.
#   Shannon, C.E. (1948). A Mathematical Theory of Communication. Bell
#     System Technical Journal, 27, 379-423.
#   Cressie, N. and Read, T.R.C. (1984). Multinomial goodness-of-fit tests.
#     Journal of the Royal Statistical Society B, 46, 440-464.

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

#' Validate a pair of compatible probability vectors
#' @keywords internal
.check_prob_pair <- function(p, q) {
  .check_prob(p, name = "p")
  .check_prob(q, name = "q")
  if (length(p) != length(q))
    stop("p and q must have the same length.")
  if (any(p > 0 & q == 0))
    stop("q must be positive wherever p is positive (support of p must be ",
         "contained in support of q).")
  invisible(NULL)
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

# ---- Cross-entropy (KL divergence) ------------------------------------------

#' Kullback-Leibler Cross-Entropy Divergence
#'
#' Computes the KL divergence D(p || q) = sum(p * log(p / q)), measuring the
#' information discrepancy when distribution q is used in place of the true
#' distribution p. The convention 0 * log(0/q) = 0 is applied.
#'
#' @param p Numeric vector. The "true" or estimated distribution.
#' @param q Numeric vector. The reference or prior distribution.  Must be
#'   positive wherever \code{p} is positive.
#' @param base Positive numeric. Logarithm base (default: natural log).
#'
#' @return A single non-negative numeric value. Returns 0 when \code{p == q}.
#'
#' @details
#' The KL divergence is not symmetric (D(p||q) != D(q||p) in general) and is
#' not a true metric. It equals zero if and only if p = q almost everywhere.
#' Under the CE (GCE) formulation, the objective is to minimize D(p||q)
#' subject to moment constraints; this is equivalent to the ME formulation
#' when q is uniform. See Golan (2008, Section 3.1, Eq. 3.4).
#'
#' @references
#' Kullback, S. and Leibler, R.A. (1951). On information and sufficiency.
#' \emph{Annals of Mathematical Statistics}, \strong{22}, 79-86.
#'
#' Golan, A. (2008). Information and Entropy Econometrics. \emph{Foundations
#' and Trends in Econometrics}, \strong{2}(1-2), 1-145.
#'
#' @examples
#' p <- c(0.5, 0.3, 0.2)
#' q <- rep(1/3, 3)          # uniform prior
#' kl_divergence(p, q)       # D(p || q) > 0
#' kl_divergence(q, q)       # = 0
#'
#' @export
kl_divergence <- function(p, q, base = exp(1)) {
  .check_prob_pair(p, q)
  if (!is.numeric(base) || length(base) != 1 || base <= 0)
    stop("base must be a single positive number.")
  idx <- p > 0
  d <- sum(p[idx] * log(p[idx] / q[idx]))
  if (base != exp(1)) d <- d / log(base)
  d
}

# ---- Renyi entropy ----------------------------------------------------------

#' Renyi Entropy of Order Alpha
#'
#' Computes the Renyi entropy of order alpha:
#' H^R_alpha(p) = (1 / (1 - alpha)) * log(sum(p^alpha)).
#' Shannon entropy is recovered in the limit as alpha -> 1.
#'
#' @param p Numeric vector of probabilities.
#' @param alpha Numeric scalar. Order parameter. Must be positive and not
#'   equal to 1. For alpha -> 1, Shannon entropy results (use
#'   \code{shannon_entropy()} directly).
#' @param base Positive numeric. Logarithm base (default: natural log).
#'
#' @return A single numeric value.
#'
#' @details
#' The Renyi entropy is a generalization of Shannon entropy indexed by a
#' single parameter alpha (Golan, 2008, Eq. 3.5). For alpha > 1, higher
#' probability events contribute more weight to the sum than they do under
#' Shannon entropy. The Renyi and Tsallis entropy measures of order (alpha+1)
#' are related to the Cressie-Read measure of order alpha via Eq. (3.9) of
#' Golan (2008).
#'
#' @references
#' Renyi, A. (1961). On measures of entropy and information. In
#' \emph{Proceedings of the 4th Berkeley Symposium}, Vol. I, 547-561.
#'
#' Golan, A. (2008). Information and Entropy Econometrics. \emph{Foundations
#' and Trends in Econometrics}, \strong{2}(1-2), 1-145.
#'
#' @examples
#' p <- c(0.5, 0.3, 0.2)
#' renyi_entropy(p, alpha = 2)
#' renyi_entropy(p, alpha = 0.5)
#'
#' # As alpha -> 1, result converges to Shannon entropy
#' renyi_entropy(p, alpha = 1 - 1e-9)
#' shannon_entropy(p)
#'
#' @export
renyi_entropy <- function(p, alpha, base = exp(1)) {
  .check_prob(p)
  if (!is.numeric(alpha) || length(alpha) != 1)
    stop("alpha must be a single numeric value.")
  if (alpha <= 0)
    stop("alpha must be positive.")
  if (abs(alpha - 1) < 1e-10)
    stop("alpha = 1 is not valid for renyi_entropy(); use shannon_entropy() instead.")
  if (!is.numeric(base) || length(base) != 1 || base <= 0)
    stop("base must be a single positive number.")

  h <- (1 / (1 - alpha)) * log(sum(p^alpha))
  if (base != exp(1)) h <- h / log(base)
  h
}

# ---- Renyi cross-entropy ----------------------------------------------------

#' Renyi Cross-Entropy Divergence of Order Alpha
#'
#' Computes the Renyi cross-entropy (divergence) of order alpha:
#' D^R_alpha(p || q) = (1 / (1 - alpha)) * log(sum(p^alpha * q^(1-alpha))).
#' Reduces to the KL divergence as alpha -> 1.
#'
#' @param p Numeric vector. The estimated distribution.
#' @param q Numeric vector. The reference (prior) distribution.
#' @param alpha Numeric scalar. Order parameter (positive, not equal to 1).
#' @param base Positive numeric. Logarithm base (default: natural log).
#'
#' @return A single numeric value.
#'
#' @references
#' Golan, A. (2008). Information and Entropy Econometrics. \emph{Foundations
#' and Trends in Econometrics}, \strong{2}(1-2), 1-145. Eq. (3.6).
#'
#' @examples
#' p <- c(0.5, 0.3, 0.2)
#' q <- rep(1/3, 3)
#' renyi_divergence(p, q, alpha = 2)
#'
#' @export
renyi_divergence <- function(p, q, alpha, base = exp(1)) {
  .check_prob_pair(p, q)
  if (!is.numeric(alpha) || length(alpha) != 1)
    stop("alpha must be a single numeric value.")
  if (alpha <= 0)
    stop("alpha must be positive.")
  if (abs(alpha - 1) < 1e-10)
    stop("alpha = 1 is not valid; use kl_divergence() instead.")
  if (!is.numeric(base) || length(base) != 1 || base <= 0)
    stop("base must be a single positive number.")

  idx <- p > 0
  d <- (1 / (1 - alpha)) * log(sum(p[idx]^alpha * q[idx]^(1 - alpha)))
  if (base != exp(1)) d <- d / log(base)
  d
}

# ---- Tsallis cross-entropy --------------------------------------------------

#' Tsallis Cross-Entropy Divergence of Order Alpha
#'
#' Computes the Tsallis cross-entropy divergence of order alpha:
#' D^T_alpha(p || q) = (1 / (1 - alpha)) * (sum(p^alpha * q^(1-alpha)) - 1).
#' Reduces to the KL divergence as alpha -> 1.
#'
#' @param p Numeric vector. The estimated distribution.
#' @param q Numeric vector. The reference (prior) distribution.
#' @param alpha Numeric scalar. Order parameter (positive, not equal to 1).
#'
#' @return A single numeric value.
#'
#' @details
#' The Tsallis measure (Golan, 2008, Eq. 3.7) is similar in functional form
#' to the Box-Cox transformation. It is pseudo-additive for independent
#' sub-systems, unlike Shannon entropy which is strictly additive.
#' The Renyi and Tsallis measures of order (alpha+1) are related to the
#' Cressie-Read measure of order alpha via Eq. (3.9) of Golan (2008).
#'
#' @references
#' Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs statistics.
#' \emph{Journal of Statistical Physics}, \strong{52}, 479-487.
#'
#' Golan, A. (2008). Information and Entropy Econometrics. \emph{Foundations
#' and Trends in Econometrics}, \strong{2}(1-2), 1-145. Eq. (3.7).
#'
#' @examples
#' p <- c(0.5, 0.3, 0.2)
#' q <- rep(1/3, 3)
#' tsallis_divergence(p, q, alpha = 2)
#'
#' @export
tsallis_divergence <- function(p, q, alpha) {
  .check_prob_pair(p, q)
  if (!is.numeric(alpha) || length(alpha) != 1)
    stop("alpha must be a single numeric value.")
  if (alpha <= 0)
    stop("alpha must be positive.")
  if (abs(alpha - 1) < 1e-10)
    stop("alpha = 1 is not valid; use kl_divergence() instead.")

  idx <- p > 0
  (1 / (1 - alpha)) * (sum(p[idx]^alpha * q[idx]^(1 - alpha)) - 1)
}

# ---- Cressie-Read -----------------------------------------------------------

#' Cressie-Read Power Divergence Statistic
#'
#' Computes the Cressie-Read power divergence of order alpha:
#' D^CR_alpha(p || q) = (1 / (alpha * (1 + alpha))) * sum(p * ((p/q)^alpha - 1)).
#'
#' This measure is the unifying criterion for the class of IT estimators.
#' Special cases include:
#' \itemize{
#'   \item alpha -> 0:  Kullback-Leibler divergence D(p||q)
#'   \item alpha -> -1: Empirical Likelihood (reverse KL)
#'   \item alpha = 1:  Pearson chi-squared / 2
#'   \item alpha = -2: Neyman chi-squared / 2
#'   \item alpha = 2/3: Freeman-Tukey statistic
#' }
#'
#' @param p Numeric vector. The estimated distribution.
#' @param q Numeric vector. The reference (prior) distribution.
#' @param alpha Numeric scalar. Order parameter. Must not be 0 or -1 (use
#'   \code{kl_divergence()} for alpha = 0 and the EL criterion for alpha = -1).
#'
#' @return A single numeric value (non-negative).
#'
#' @details
#' The Cressie-Read (CR) family (Golan, 2008, Eq. 3.8) provides the
#' objective function for all IT estimators in the Golan framework. Different
#' values of alpha select different estimators within the generic IT class.
#' Only alpha -> 0 (Shannon / KL) is fully consistent with information theory.
#' The Renyi and Tsallis measures of order (alpha+1) are related to the CR
#' measure of order alpha via the unifying equation (Golan, 2008, Eq. 3.9):
#'
#'   D^R_{alpha+1}(p||q) = -(1/alpha) * log[1 + alpha*(alpha+1)*D^CR_alpha(p||q)]
#'
#' @references
#' Cressie, N. and Read, T.R.C. (1984). Multinomial goodness-of-fit tests.
#' \emph{Journal of the Royal Statistical Society B}, \strong{46}, 440-464.
#'
#' Golan, A. (2008). Information and Entropy Econometrics. \emph{Foundations
#' and Trends in Econometrics}, \strong{2}(1-2), 1-145. Eq. (3.8), (3.9).
#'
#' @examples
#' p <- c(0.5, 0.3, 0.2)
#' q <- rep(1/3, 3)
#'
#' # alpha = 1 gives (Pearson chi-squared) / 2
#' cressie_read(p, q, alpha = 1)
#'
#' # Verify connection to KL divergence (alpha -> 0)
#' cressie_read(p, q, alpha = 1e-9)
#' kl_divergence(p, q)
#'
#' @export
cressie_read <- function(p, q, alpha) {
  .check_prob_pair(p, q)
  if (!is.numeric(alpha) || length(alpha) != 1)
    stop("alpha must be a single numeric value.")
  if (abs(alpha) < 1e-10)
    stop("alpha = 0 is not valid for cressie_read(); use kl_divergence() instead.")
  if (abs(alpha + 1) < 1e-10)
    stop("alpha = -1 corresponds to the EL criterion; use el() instead.")

  idx <- p > 0
  (1 / (alpha * (1 + alpha))) * sum(p[idx] * ((p[idx] / q[idx])^alpha - 1))
}

# ---- Normalized entropy -----------------------------------------------------

#' Normalized Entropy (Signal Measure)
#'
#' Computes the normalized entropy S(p) = H(p) / H(q), where H is the Shannon
#' entropy. When q is uniform (the default), this simplifies to
#' S(p) = H(p) / log(K), where K = length(p).
#'
#' @param p Numeric vector of probabilities (the estimated distribution).
#' @param q Numeric vector. Reference (prior) distribution. Defaults to
#'   the uniform distribution over \code{length(p)} outcomes.
#'
#' @return A numeric value in [0, 1]. Returns 1 when p is uniform (complete
#'   ignorance) and 0 when p is a point mass (perfect certainty).
#'
#' @details
#' The normalized entropy S(p) is the primary goodness-of-fit statistic for
#' IT models (Golan, 2008, Section 6.4). It measures the proportion of
#' maximum uncertainty that remains after conditioning on the data. A value
#' close to 1 indicates the data carry little information; a value close to 0
#' indicates the model fits well. The entropy-ratio test statistic and
#' pseudo-R-squared are direct functions of S(p).
#'
#' For non-uniform priors q, the normalized entropy is defined as
#' S(p) = H(p) / H(q), where H(p) = -sum(p * log(p)).
#'
#' @references
#' Golan, A. (2008). Information and Entropy Econometrics. \emph{Foundations
#' and Trends in Econometrics}, \strong{2}(1-2), 1-145. Section 6.4.
#'
#' @examples
#' # Uniform distribution -> S = 1 (complete ignorance)
#' normalized_entropy(rep(1/6, 6))
#'
#' # Point mass -> S = 0 (perfect certainty)
#' normalized_entropy(c(1, 0, 0, 0, 0, 0))
#'
#' # Estimated probabilities after ME estimation
#' p_hat <- c(0.05, 0.10, 0.20, 0.35, 0.25, 0.05)
#' normalized_entropy(p_hat)
#'
#' @export
normalized_entropy <- function(p, q = NULL) {
  .check_prob(p)
  K <- length(p)
  if (is.null(q)) {
    h_max <- log(K)
  } else {
    .check_prob(q, name = "q")
    if (length(q) != K) stop("p and q must have the same length.")
    h_max <- shannon_entropy(q)
  }
  if (h_max < .Machine$double.eps)
    stop("H(q) is zero; the reference distribution is degenerate.")
  shannon_entropy(p) / h_max
}
