# utils.R
# Shared utility functions for the infometrics package.

#' Normalize a numeric vector to [0, 1]
#'
#' Rescales a numeric vector so that all values lie in [0, 1] by dividing by
#' the maximum absolute value. Useful for pre-scaling moment matrices before
#' ME/CE estimation to avoid numerical overflow in the exponential terms of
#' the partition function.
#'
#' @param x Numeric vector or matrix.
#' @param by Character. Either \code{"max"} (divide by \code{max(abs(x))},
#'   default) or \code{"range"} (min-max scaling to [0, 1]).
#'
#' @return A numeric vector or matrix of the same dimensions as \code{x},
#'   scaled to [0, 1].
#'
#' @details
#' Golan (2008, Section 7.1) notes that normalization is often necessary when
#' data contain exponential terms (as in the partition function Omega), since
#' the concentrated (dual) ME/CE model involves expressions of the form
#' exp(lambda * x). A common normalization is dividing each element by
#' max\{x_j, y_i\}.
#'
#' @examples
#' x <- c(10, 30, 50, 20, 40)
#' normalize_data(x)
#' normalize_data(x, by = "range")
#'
#' @export
normalize_data <- function(x, by = c("max", "range")) {
  by <- match.arg(by)
  if (!is.numeric(x)) stop("x must be numeric.")
  if (by == "max") {
    m <- max(abs(x))
    if (m < .Machine$double.eps) return(x)
    x / m
  } else {
    rng <- range(x)
    if (diff(rng) < .Machine$double.eps) return(x - rng[1])
    (x - rng[1]) / diff(rng)
  }
}

#' Construct a Symmetric Support Space
#'
#' Creates an M-point symmetric support vector centered at zero (or a
#' specified center) over a given half-range. Support spaces are required
#' inputs for the GME and GCE estimators, bounding the parameter vector beta
#' and the error vector epsilon.
#'
#' @param half_range Positive numeric. The half-width of the support interval.
#'   The support spans from \code{center - half_range} to
#'   \code{center + half_range}.
#' @param M Positive integer. Number of support points. Must be >= 2.
#'   The default of 5 (a 5-point support) is commonly used in the GME
#'   literature.
#' @param center Numeric. Center of the support interval (default 0).
#'
#' @return A numeric vector of length M giving the support points, evenly
#'   spaced from \code{center - half_range} to \code{center + half_range}.
#'
#' @details
#' In the GME framework (Golan, 2008, Section 6.1), each parameter beta_k is
#' reparameterized as the expected value of a random variable defined on a
#' bounded support z_k = (z_\{k1\}, ..., z_\{kM\}). The support must contain the
#' true parameter value. A common default is a symmetric 5-point grid with
#' half-range set to 3 standard deviations of the OLS residuals (for the error
#' support V) or 2 times the absolute OLS estimate (for the signal support Z).
#'
#' The complexity of the GME dual model is invariant to M (the number of
#' support points), since the real parameters are the Lagrange multipliers
#' whose dimension equals T (the number of observations), not K*M.
#'
#' @references
#' Golan, A. (2008). Information and Entropy Econometrics. \emph{Foundations
#' and Trends in Econometrics}, \strong{2}(1-2), 1-145. Section 6.1.
#'
#' @examples
#' # 5-point support for error terms: [-3, -1.5, 0, 1.5, 3]
#' make_support(half_range = 3, M = 5)
#'
#' # 3-point support centered at a non-zero value: [0, 5, 10]
#' make_support(half_range = 5, M = 3, center = 5)
#'
#' # 7-point support for a coefficient bounded in [-2, 2]
#' make_support(half_range = 2, M = 7)
#'
#' @export
make_support <- function(half_range, M = 5L, center = 0) {
  if (!is.numeric(half_range) || length(half_range) != 1 || half_range <= 0)
    stop("half_range must be a single positive number.")
  M <- as.integer(M)
  if (M < 2L)
    stop("M must be at least 2.")
  if (!is.numeric(center) || length(center) != 1)
    stop("center must be a single numeric value.")
  seq(from = center - half_range, to = center + half_range, length.out = M)
}

#' Default Support Spaces from OLS
#'
#' Constructs default signal (Z) and error (V) support spaces for GME
#' estimation using OLS estimates as reference points. The signal support for
#' each coefficient is centered on the OLS estimate with half-range equal to
#' \code{signal_scale} times the absolute OLS estimate. The error support is
#' symmetric around zero with half-range equal to \code{error_scale} times
#' the OLS residual standard deviation.
#'
#' @param y Numeric vector of length T. The dependent variable.
#' @param X Numeric matrix of dimension T x K. The regressor matrix.
#' @param M_signal Integer. Number of support points per coefficient (default 5).
#' @param M_error Integer. Number of support points for error terms (default 5).
#' @param signal_scale Positive numeric. Half-range multiplier for the signal
#'   support (default 2: support spans +/- 2 * |OLS estimate|).
#' @param error_scale Positive numeric. Half-range multiplier for the error
#'   support (default 3: support spans +/- 3 * OLS residual SD, the
#'   "three-sigma rule" of Pukelsheim, 1994).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{Z}}{K x M_signal matrix. Row k is the support for beta_k.}
#'     \item{\code{V}}{Numeric vector of length M_error. The error support,
#'       symmetric around zero.}
#'   }
#'
#' @references
#' Golan, A. (2008). Information and Entropy Econometrics. \emph{Foundations
#' and Trends in Econometrics}, \strong{2}(1-2), 1-145. Section 6.1.
#'
#' Pukelsheim, F. (1994). The three sigma rule.
#' \emph{The American Statistician}, \strong{48}, 88-91.
#'
#' @examples
#' set.seed(42)
#' n <- 50; K <- 3
#' X <- cbind(1, matrix(rnorm(n * (K-1)), n, K-1))
#' y <- X %*% c(1, 2, -1) + rnorm(n)
#' sp <- default_supports(y, X)
#' sp$Z   # 3 x 5 signal support matrix
#' sp$V   # 5-point error support
#'
#' @export
default_supports <- function(y, X, M_signal = 5L, M_error = 5L,
                             signal_scale = 2, error_scale = 3) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  K <- ncol(X)
  T_obs <- nrow(X)

  if (length(y) != T_obs)
    stop("length(y) must equal nrow(X).")
  if (M_signal < 2L) stop("M_signal must be >= 2.")
  if (M_error  < 2L) stop("M_error must be >= 2.")

  # OLS estimates
  ols <- tryCatch(
    as.numeric(solve(crossprod(X), crossprod(X, y))),
    error = function(e) {
      warning("OLS failed (singular X'X); using zero-centered signal supports.")
      rep(0, K)
    }
  )
  ols_resid_sd <- stats::sd(y - as.numeric(X %*% ols))
  if (!is.finite(ols_resid_sd) || ols_resid_sd < .Machine$double.eps)
    ols_resid_sd <- 1.0

  # Signal support: centered on OLS estimate, half-range = scale * |beta_ols|
  Z <- matrix(0, nrow = K, ncol = M_signal)
  for (k in seq_len(K)) {
    hr <- max(signal_scale * abs(ols[k]), 1e-4)  # minimum half-range
    Z[k, ] <- make_support(half_range = hr, M = M_signal, center = ols[k])
  }

  # Error support: symmetric around zero, half-range = scale * sigma_ols
  V <- make_support(half_range = error_scale * ols_resid_sd, M = M_error)

  list(Z = Z, V = V)
}

# ---- Fano error bounds for a row-stochastic probability matrix --------------
# Model-agnostic (Golan 2008, sec 3.6/7.5): each row of p is a distribution over
# K categories; the modal classifier's error pe = 1 - max_k p_k is lower-bounded
# by S(p) - log(2)/log(K), where S(p) = H(p)/log(K) is the normalized entropy.
# Returns a per-row data frame with an "overall" attribute.

#' @keywords internal
.fano_row_bounds <- function(p, K) {
  p_max <- apply(p, 1L, max)
  pe    <- 1 - p_max
  H     <- -rowSums(ifelse(p > 0, p * log(p), 0))          # nats
  S     <- H / log(K)
  pe_lo <- pmax(0, S - log(2) / log(K))                    # Golan sec 7.5 weak bound
  per   <- data.frame(p_max = p_max, pe = pe, H = H, S = S, pe_lower = pe_lo)
  attr(per, "overall") <- c(mean_pe = mean(pe), mean_pe_lower = mean(pe_lo),
                            S_system = sum(H) / (nrow(p) * log(K)))
  per
}
