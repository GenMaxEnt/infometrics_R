# linreg_iv.R
# (Relaxed) stochastic-moments GME estimator for instrumental-variables
# regression (Golan 2008, pp. 89-91). IV sibling of linreg().
#
# Model y = X beta + e, identified by instrument moments IV'(y - X beta - e) = 0,
# with beta = Z p (signal support Z, K x M) and e = v w (error support v),
# weighted by nu in (0, 1). The concentrated dual is solved over the Lagrange
# multipliers lambda (one per instrument moment):
#
#   p_km prop to p0_km exp( z_km (X'IV lambda)_k / nu )          (signal)
#   w_ij prop to w0_ij exp( (IV lambda)_i v_j / (1 - nu) )       (noise)
#   e_i  = sum_j v_j w_ij ;  beta_k = sum_m z_km p_km
#   dual (MAXIMISED): y'IV lambda - nu sum_k logOmega_k - (1-nu) sum_i logPsi_i
#   gradient:          IV'(y - X beta - e)
#
# Notes:
#  - The dual is kept in the maximise form (fnscale = -1) by user choice (with
#    markov_ce/markov_gce/multinomial_gce).
#  - The instruments are standardized internally for conditioning (beta is
#    invariant to instrument scaling); lambda is reported on the original scale.
#  - Supports just-identified (ncol(IV) == ncol(X)) and over-identified
#    (ncol(IV) > ncol(X)) systems.
#  - Standard errors are deferred (the GME-IV asymptotic variance is non-trivial).
#
# References:
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145. Pp. 89-91.

# ---- internal: stable signal / noise pieces --------------------------------

#' @keywords internal
.linreg_iv_pieces <- function(lambda, XtIV, IV, Z, logp0, v, logw0, nu) {
  a <- as.vector(XtIV %*% lambda)                 # K  = (X'IV) lambda
  L <- logp0 + Z * (a / nu)                       # K x M : z_km a_k / nu
  rmax <- apply(L, 1L, max); EL <- exp(L - rmax); rs <- rowSums(EL)
  p <- EL / rs; b <- rowSums(Z * p)               # signal probs + beta
  g <- as.vector(IV %*% lambda)                   # N  = (IV lambda)_i
  B <- logw0 + outer(g, v) / (1 - nu)             # N x J
  cmax <- apply(B, 1L, max); EB <- exp(B - cmax); cs <- rowSums(EB)
  w <- EB / cs; e <- as.vector(w %*% v)           # noise weights + e
  list(p = p, b = b, w = w, e = e,
       logOmega = rmax + log(rs), logPsi = cmax + log(cs))
}

#' @keywords internal
.linreg_iv_obj <- function(lambda, y, X, IV, XtIV, IVty, Z, logp0, v, logw0, nu) {
  pc <- .linreg_iv_pieces(lambda, XtIV, IV, Z, logp0, v, logw0, nu)
  sum(IVty * lambda) - nu * sum(pc$logOmega) - (1 - nu) * sum(pc$logPsi)
}

#' @keywords internal
.linreg_iv_grad <- function(lambda, y, X, IV, XtIV, IVty, Z, logp0, v, logw0, nu) {
  pc <- .linreg_iv_pieces(lambda, XtIV, IV, Z, logp0, v, logw0, nu)
  as.vector(crossprod(IV, y - X %*% pc$b - pc$e))
}

# ===========================================================================
#  linreg_iv()
# ===========================================================================

#' Stochastic-Moments GME Estimator for IV Regression
#'
#' @description
#' Fits the linear model \eqn{y = X\beta + e} by the (relaxed) stochastic-moments
#' Generalized Maximum Entropy estimator of Golan (2008, pp. 89-91), identified
#' by instrument moments \eqn{IV'(y - X\beta - e) = 0}. Each coefficient is
#' reparameterized on a bounded signal support, \eqn{\beta_k = \sum_m z_{km}
#' p_{km}}, and each error on a noise support, \eqn{e_i = \sum_j v_j w_{ij}}; the
#' signal and noise entropies are maximized subject to the instrument moments,
#' with weight \eqn{\nu \in (0,1)}. It is the instrumental-variables sibling of
#' \code{\link{linreg}}.
#'
#' @details
#' The concentrated dual is maximized over the Lagrange multipliers \eqn{\lambda}
#' (one per instrument moment): with \eqn{p_{km} \propto p0_{km}\exp(z_{km}
#' (X'IV\lambda)_k/\nu)} and \eqn{w_{ij} \propto w0_{ij}\exp((IV\lambda)_i
#' v_j/(1-\nu))}, the gradient is the instrument moment
#' \eqn{IV'(y - X\beta - e)}. Supports just-identified
#' (\code{ncol(IV) == ncol(X)}) and over-identified (\code{ncol(IV) > ncol(X)})
#' systems. The instruments are standardized internally for numerical
#' conditioning (\eqn{\beta} is invariant to instrument scaling); \code{lambda}
#' is reported on the original scale. \code{coef()} returns \eqn{\beta}.
#'
#' As the signal support widens the GME estimate approaches the exact (2SLS-type)
#' IV solution; narrower supports shrink \eqn{\beta} toward the support centers.
#' Standard errors are not currently reported.
#'
#' @param y Numeric response vector of length N.
#' @param X N-by-K design matrix (include an intercept column if wanted).
#' @param IV N-by-P instrument matrix, \code{P >= K}. For exogenous regressors,
#'   use the regressor as its own instrument.
#' @param Z Coefficient (signal) support: a K-by-M matrix whose row k is the
#'   support for \eqn{\beta_k} (as in \code{\link{linreg}}).
#' @param p0 Optional K-by-M signal prior (rows strictly positive, summing to 1);
#'   default uniform.
#' @param v Optional error support: \code{NULL} (default 3-point symmetric grid
#'   on \eqn{\pm 3\,\mathrm{sd}(y)}), a single whole number giving the number of
#'   support points, or an explicit numeric vector.
#' @param w0 Optional N-by-J error prior (rows strictly positive, summing to 1);
#'   default uniform.
#' @param nu Signal/noise entropy weight in \eqn{(0,1)}; default 0.5.
#' @param control Named list merged over defaults and passed to
#'   \code{\link[stats]{optim}} (BFGS): \code{maxit} (default 1000) and
#'   \code{reltol} (default 1e-12). \code{fnscale} is forced to -1 (the dual is
#'   maximized).
#'
#' @return An object of class \code{c("linreg_iv", "infometrics")} with
#'   \code{coefficients}/\code{b_hat} (\eqn{\beta}, length K; \code{coef()}
#'   returns it), \code{lambda_hat} (length P, original instrument scale),
#'   \code{p_hat} (K-by-M), \code{w_hat} (N-by-J), \code{e_hat}/\code{e} (N),
#'   \code{fitted.values} (\eqn{X\beta}), \code{residuals} (\eqn{y - X\beta}),
#'   \code{H_signal} (K-vector), \code{S}, \code{objective}/\code{value},
#'   \code{moment_resid}, \code{converged}/\code{convergence}, \code{method},
#'   plus \code{Z}, \code{p0}, \code{v}, \code{w0}, \code{nu}, \code{X},
#'   \code{y}, \code{IV}, \code{N}/\code{K}/\code{P}/\code{M}, \code{call}.
#'
#' @references Golan, A. (2008). \emph{Information and Entropy Econometrics -
#'   A Review and Synthesis}. Foundations and Trends in Econometrics, 2(1-2),
#'   1-145. Pages 89-91.
#'
#' @seealso \code{\link{linreg}} for the (non-IV) GME/GCE regression.
#'
#' @examples
#' set.seed(1)
#' n  <- 200L
#' z  <- rnorm(n)                         # instrument
#' u  <- rnorm(n)                         # structural error
#' xe <- 0.7 * z + u + rnorm(n)           # endogenous regressor (corr. with u)
#' y  <- 1 + 1.5 * xe + u                 # true slope 1.5; OLS biased
#' X  <- cbind(1, xe)                     # intercept + endogenous regressor
#' IV <- cbind(1, z)                      # intercept is its own instrument
#' Z  <- matrix(c(-10, 0, 10), nrow = 2, ncol = 3, byrow = TRUE)  # signal support
#' fit <- linreg_iv(y, X, IV, Z)
#' coef(fit)                              # slope ~ 1.5 (vs OLS biased upward)
#'
#' @importFrom stats optim sd
#' @export
linreg_iv <- function(y, X, IV, Z, p0 = NULL, v = NULL, w0 = NULL, nu = 0.5,
                      control = list()) {
  y <- as.numeric(y); X <- as.matrix(X); IV <- as.matrix(IV); Z <- as.matrix(Z)
  N <- nrow(X); K <- ncol(X); P <- ncol(IV); M <- ncol(Z)

  if (length(y) != N) stop("length(y) (", length(y), ") must equal nrow(X) (", N, ").")
  if (nrow(IV) != N) stop("nrow(IV) (", nrow(IV), ") must equal nrow(X) (", N, ").")
  if (P < K) stop("Under-identified: ncol(IV) (", P, ") must be >= ncol(X) (", K, ").")
  if (nrow(Z) != K) stop("nrow(Z) (", nrow(Z), ") must equal ncol(X) (", K, ").")
  if (!is.numeric(nu) || length(nu) != 1L || nu <= 0 || nu >= 1)
    stop("'nu' must be a single number in (0, 1). Got nu = ", nu, ".")

  if (is.null(p0)) p0 <- matrix(1 / M, K, M)
  else {
    p0 <- as.matrix(p0)
    if (!all(dim(p0) == c(K, M))) stop("dim(p0) must match dim(Z) = ", K, " x ", M, ".")
    if (any(p0 <= 0)) stop("All signal prior masses 'p0' must be strictly positive.")
    if (any(abs(rowSums(p0) - 1) > 1e-8)) stop("Rows of 'p0' must sum to 1.")
  }

  sdY <- stats::sd(y)
  if (is.null(v)) {
    j <- if (is.null(w0)) 3L else ncol(w0)
    v <- seq(-3 * sdY, 3 * sdY, length.out = j)
  } else if (length(v) == 1L) {
    j <- as.integer(v); if (j < 2L) stop("scalar 'v' (support count) must be >= 2.")
    v <- seq(-3 * sdY, 3 * sdY, length.out = j)
  } else { v <- as.numeric(v); j <- length(v) }

  if (is.null(w0)) w0 <- matrix(1 / j, N, j)
  else {
    w0 <- as.matrix(w0)
    if (!all(dim(w0) == c(N, j))) stop("dim(w0) must be ", N, " x ", j, ".")
    if (any(w0 <= 0)) stop("All noise prior masses 'w0' must be strictly positive.")
    if (any(abs(rowSums(w0) - 1) > 1e-8)) stop("Rows of 'w0' must sum to 1.")
  }

  # ---- standardize instruments for conditioning (beta invariant) -----------
  iv_scale <- apply(IV, 2L, function(z) { s <- stats::sd(z); if (s > 0) s else 1 })
  IVs <- sweep(IV, 2L, iv_scale, "/")

  logp0 <- log(p0); logw0 <- log(w0)
  XtIV <- crossprod(X, IVs)            # K x P  (constant)
  IVty <- as.vector(crossprod(IVs, y)) # P      (constant)

  con <- list(maxit = 1000L, reltol = 1e-12)
  con[names(control)] <- control
  con$fnscale <- -1

  est <- stats::optim(
    par = rep(0, P), fn = .linreg_iv_obj, gr = .linreg_iv_grad,
    y = y, X = X, IV = IVs, XtIV = XtIV, IVty = IVty,
    Z = Z, logp0 = logp0, v = v, logw0 = logw0, nu = nu,
    method = "BFGS", control = con
  )
  if (est$convergence != 0L)
    warning("optim did not converge (code ", est$convergence, ").")

  pc <- .linreg_iv_pieces(est$par, XtIV, IVs, Z, logp0, v, logw0, nu)
  b  <- pc$b
  lambda <- est$par / iv_scale          # report lambda on the original IV scale

  rn <- colnames(X); if (is.null(rn)) rn <- paste0("x", seq_len(K))
  names(b) <- rn; rownames(pc$p) <- rn
  yhat <- as.vector(X %*% b); resid <- y - yhat

  ent <- function(z) { z <- z[z > 0]; -sum(z * log(z)) }
  H_signal <- apply(pc$p, 1L, ent)
  S <- mean(H_signal) / log(M)
  moment_resid <- max(abs(crossprod(IVs, y - X %*% b - pc$e)))

  structure(
    list(
      # draft names
      lambda_hat = lambda, p_hat = pc$p, b_hat = b, w_hat = pc$w, e_hat = pc$e,
      convergence = est$convergence,
      # canonical aliases
      coefficients = b, fitted.values = yhat, residuals = resid,
      e = pc$e, H_signal = H_signal, S = S,
      objective = est$value, value = est$value,
      converged = (est$convergence == 0L), method = "dual",
      # extras
      Z = Z, p0 = p0, v = v, w0 = w0, nu = nu, X = X, y = y, IV = IV,
      iv_scale = iv_scale, moment_resid = moment_resid,
      N = N, K = K, P = P, M = M, call = match.call()
    ),
    class = c("linreg_iv", "infometrics")
  )
}

# ===========================================================================
#  S3 methods
# ===========================================================================

#' @export
coef.linreg_iv <- function(object, ...) object$coefficients

#' @export
fitted.linreg_iv <- function(object, ...) object$fitted.values

#' @export
residuals.linreg_iv <- function(object, ...) object$residuals

#' @export
print.linreg_iv <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  ident <- if (x$P == x$K) "just-identified" else "over-identified"
  cat("Stochastic-moments GME-IV (Golan 2008, pp. 89-91)\n")
  cat(strrep("-", 50L), "\n")
  cat(sprintf("  N = %d   K = %d   instruments = %d (%s)   M = %d   nu = %.3g\n",
              x$N, x$K, x$P, ident, x$M, x$nu))
  cat(sprintf("  signal support : [%.4g, %.4g]   error support : [%.4g, %.4g]\n",
              min(x$Z), max(x$Z), min(x$v), max(x$v)))
  cat(sprintf("  normalized S = %.4f   max|moment resid| = %s   convergence = %s\n",
              x$S, format(x$moment_resid, digits = digits),
              if (x$converged) "yes" else "no"))
  cat("\nCoefficients (beta):\n")
  print(round(x$coefficients, digits))
  invisible(x)
}

#' @export
summary.linreg_iv <- function(object, digits = max(3L, getOption("digits") - 3L),
                              ...) {
  ident <- if (object$P == object$K) "just-identified" else "over-identified"
  cat("\nStochastic-moments GME-IV (Golan 2008, pp. 89-91)\n")
  cat("Call: "); print(object$call); cat("\n")
  cat(sprintf("  N=%d  K=%d  instruments=%d (%s)  M=%d  nu=%.3g\n",
              object$N, object$K, object$P, ident, object$M, object$nu))
  cat("\nResiduals (y - X beta):\n")
  print(summary(object$residuals, digits = digits))
  cat("\nCoefficients (beta = Z p):\n")
  print(round(object$coefficients, digits))
  cat(sprintf("\nnormalized signal entropy S = %.4f   max|moment resid| = %s\n",
              object$S, format(object$moment_resid, digits = digits)))
  cat(sprintf("convergence: %s\n", if (object$converged) "yes (0)" else "no"))
  cat("\nNote: standard errors are not reported; the GME-IV asymptotic variance\n",
      "is deferred to the package inference layer.\n", sep = "")
  invisible(object)
}
