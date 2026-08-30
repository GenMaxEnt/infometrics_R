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

# ---- standard-error internals ----------------------------------------------
# beta = Z p(lambda) is a Z-estimator of the instrument moments IV'(y-Xb-e)=0,
# so Var(beta) = Jb A^{-1} Omega A^{-1} Jb', with A the dual Hessian,
# Jb = dbeta/dlambda, and Omega the "meat". All on the standardized IV scale
# (beta is scale-invariant). Validated vs Monte Carlo and the 2SLS robust-SE
# limit.

#' @keywords internal
.linreg_iv_analytic_vcov <- function(X, IVs, Z, v, nu, p, w, beta, e, r, meat) {
  K <- ncol(X); Tn <- nrow(X)
  varp <- rowSums(Z^2 * p) - beta^2                 # K : Var_p(z_k)
  varw <- as.vector(w %*% (v^2)) - e^2              # N : Var_w(v)_i
  XtIVs <- crossprod(X, IVs)                        # K x P
  Jb <- (varp / nu) * XtIVs                         # K x P = dbeta / dlambda
  A  <- t(XtIVs) %*% ((varp / nu) * XtIVs) +
        crossprod(IVs, (varw / (1 - nu)) * IVs)     # P x P dual Hessian (PD)
  Ai <- solve(A)
  Omega <- if (identical(meat, "sandwich")) crossprod(IVs, (r^2) * IVs)
           else (sum(r^2) / max(Tn - K, 1L)) * crossprod(IVs)
  V <- Jb %*% (Ai %*% Omega %*% Ai) %*% t(Jb)       # K x K
  V <- (V + t(V)) / 2                               # symmetrize
  dimnames(V) <- list(names(beta), names(beta))
  V
}

#' @keywords internal
.linreg_iv_boot_vcov <- function(object, B) {
  y <- object$y; X <- object$X; IV <- object$IV; n <- length(y)
  bs <- matrix(NA_real_, B, object$K)
  for (bb in seq_len(B)) {
    id <- sample.int(n, n, replace = TRUE)
    f <- tryCatch(
      linreg_iv(y[id], X[id, , drop = FALSE], IV[id, , drop = FALSE],
                Z = object$Z, p0 = object$p0, v = object$v, nu = object$nu,
                se_method = "none"),
      error = function(e) NULL)
    if (!is.null(f)) bs[bb, ] <- f$coefficients
  }
  V <- stats::cov(bs, use = "complete.obs")
  dimnames(V) <- list(names(object$coefficients), names(object$coefficients))
  V
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
#'
#' \strong{Standard errors.} Because \eqn{\hat\lambda} solves the instrument
#' moments, \eqn{\hat\beta = Zp(\hat\lambda)} is a Z-estimator with sandwich
#' covariance \eqn{\mathrm{Var}(\hat\beta) = J_\beta A^{-1}\Omega A^{-1}
#' J_\beta'}, where \eqn{A} is the dual Hessian and \eqn{J_\beta =
#' \partial\beta/\partial\lambda}. The "meat" \eqn{\Omega} is set by
#' \code{se_method}: \code{"sandwich"} (default, robust HC0
#' \eqn{IV'\mathrm{diag}(r^2)IV}), \code{"delta"} (classical homoskedastic
#' \eqn{\hat\sigma^2 IV'IV}), or \code{"bootstrap"} (pairs resampling). The
#' robust sandwich matches a Monte-Carlo sampling SD and reduces to the 2SLS
#' robust SE as the support widens; the classical delta assumes homoskedasticity.
#' These are \emph{asymptotic} SEs: the sampling distribution of \eqn{\hat\beta}
#' is support-bounded and can be skewed, so symmetric Wald intervals are
#' approximate with weak instruments or small \eqn{n}.
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
#'   on \eqn{\pm 3\,\mathrm{sd}(y)}, the three-sigma rule), a single whole number
#'   giving the number of support points, or an explicit numeric vector. Unlike
#'   \code{\link{linreg}}, this default does \emph{not} need to widen with the
#'   sample size; see \emph{Choosing the supports}.
#' @param w0 Optional N-by-J error prior (rows strictly positive, summing to 1);
#'   default uniform.
#' @param nu Entropy weight in \eqn{(0,1)} on the noise term, hence the weight
#'   on the signal is \eqn{1-\nu}; default 0.5 (equal weight on both the signal
#'   and the noise).
#' @param control Named list merged over defaults and passed to
#'   \code{\link[stats]{optim}} (BFGS): \code{maxit} (default 1000) and
#'   \code{reltol} (default 1e-12). \code{fnscale} is forced to -1 (the dual is
#'   maximized).
#' @param se_method Standard-error method: \code{"sandwich"} (default, robust),
#'   \code{"delta"} (classical), \code{"bootstrap"} (pairs resampling), or
#'   \code{"none"} (skip; \code{se_beta}/\code{vcov} left \code{NULL}).
#' @param boot Number of bootstrap resamples when
#'   \code{se_method = "bootstrap"} (default 200).
#'
#' @return
#' An object of class \code{c("linreg_iv", "infometrics")}, which is a list
#' containing the following components:
#' \describe{
#'   \item{\code{coefficients}, \code{b_hat}}{Numeric vector of length K: the
#'     estimated coefficients \eqn{\hat\beta = Z\hat p}{beta-hat = Z p-hat}
#'     (two names for the same object). Extracted by \code{\link{coef}}.}
#'   \item{\code{lambda_hat}}{Numeric vector of length P: the estimated Lagrange
#'     multipliers, one per instrument moment, reported on the original
#'     instrument scale.}
#'   \item{\code{p_hat}}{K-by-M matrix of estimated signal probabilities.}
#'   \item{\code{w_hat}}{N-by-J matrix of estimated noise probabilities.}
#'   \item{\code{e_hat}, \code{e}}{Numeric vector of length N: the estimated
#'     noise, \eqn{e_i = \sum_j v_j w_{ij}} (two names for the same object).}
#'   \item{\code{fitted.values}}{Numeric vector of length N:
#'     \eqn{X\hat\beta}{X beta-hat}. Extracted by \code{\link{fitted}}.}
#'   \item{\code{residuals}}{Numeric vector of length N:
#'     \eqn{y - X\hat\beta}{y - X beta-hat}. Extracted by
#'     \code{\link{residuals}}.}
#'   \item{\code{se_beta}}{Numeric vector of length K: standard errors of
#'     \eqn{\hat\beta}{beta-hat}, or \code{NULL} when
#'     \code{se_method = "none"}.}
#'   \item{\code{vcov}}{K-by-K covariance matrix of \eqn{\hat\beta}{beta-hat},
#'     or \code{NULL} when \code{se_method = "none"}. Extracted by
#'     \code{\link{vcov}}.}
#'   \item{\code{se_method}, \code{boot}}{The standard-error method used and,
#'     for the bootstrap, the number of resamples.}
#'   \item{\code{H_signal}}{Numeric vector of length K: the per-coefficient
#'     signal entropies.}
#'   \item{\code{S}}{Single number: the normalized signal entropy.}
#'   \item{\code{objective}, \code{value}}{Single number: the maximized dual
#'     objective (two names for the same value).}
#'   \item{\code{moment_resid}}{Single number: the largest absolute instrument
#'     moment at the optimum. Not normalized by \eqn{N} -- see
#'     \emph{Note on moment_resid}.}
#'   \item{\code{converged}}{Logical: \code{TRUE} when \code{optim} reported
#'     convergence.}
#'   \item{\code{convergence}}{Integer: the raw \code{\link[stats]{optim}}
#'     convergence code.}
#'   \item{\code{method}}{Character string naming the solver, \code{"dual"}.}
#'   \item{\code{Z}, \code{p0}, \code{v}, \code{w0}, \code{nu}, \code{X},
#'     \code{y}, \code{IV}}{The resolved inputs.}
#'   \item{\code{N}, \code{K}, \code{P}, \code{M}}{Integers: the numbers of
#'     observations, coefficients, instruments and signal support points.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @section Choosing the supports:
#' The \strong{signal} support \code{Z} is the consequential input here. It must
#' contain the true coefficients: if it does not, \eqn{\hat\beta} is pinned to
#' the boundary of \code{Z} (for example, a support of \eqn{\pm 1} returns
#' \eqn{\hat\beta_k = 1} for a coefficient whose true value is 1.5). Widening
#' \code{Z} moves the estimate toward the exact (2SLS-type) IV solution;
#' narrowing it shrinks \eqn{\beta} toward the support centers.
#'
#' The \strong{error} support \code{v} is far less critical, and -- unlike
#' \code{\link{linreg}} and \code{\link{inverse_noise}} -- its default does not
#' need to grow with the sample size. Those estimators impose one dual condition
#' per observation, \eqn{y_t = x_t'\beta + e_t}, so every realized error must fit
#' inside \code{v}; because the largest of \eqn{T} errors grows like
#' \eqn{\sigma\sqrt{2\log T}}, a fixed three-sigma support eventually becomes
#' infeasible and their duals become unbounded. \code{linreg_iv()} instead
#' imposes only the \eqn{P} aggregate instrument moments
#' \eqn{IV'(y - X\beta - e) = 0}, so the noise never has to absorb individual
#' residuals and no per-observation feasibility condition arises. In practice
#' \eqn{\hat\beta} is nearly invariant to the width of \code{v}, and large
#' samples pose no difficulty: the multipliers stay small and the estimate
#' tracks the exact IV solution.
#'
#' @section Note on \code{moment_resid}:
#' \code{moment_resid} is \eqn{\max_p |IV'(y - X\beta - e)|_p} evaluated on the
#' standardized instruments. It is a \emph{sum} over the \eqn{N} observations and
#' is not divided by \eqn{N}, so its magnitude grows with the sample size even
#' when the fit is excellent; judge it relative to \eqn{N} (or compare fits of
#' the same size) rather than against a fixed threshold.
#'
#' @references Golan, A. (2008). \emph{Information and Entropy Econometrics -
#'   A Review and Synthesis}. Foundations and Trends in Econometrics, 2(1-2),
#'   1-145. Pages 89-91. Golan, A., Judge, G. and Miller, D. (1996).
#'   \emph{Maximum Entropy Econometrics: Robust Estimation with Limited Data}.
#'   Wiley.
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
#' @importFrom stats optim sd cov pnorm printCoefmat vcov
#' @export
linreg_iv <- function(y, X, IV, Z, p0 = NULL, v = NULL, w0 = NULL, nu = 0.5,
                      se_method = c("sandwich", "delta", "bootstrap", "none"),
                      control = list(), boot = 200L) {
  se_method <- match.arg(se_method); boot <- as.integer(boot)
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

  out <- structure(
    list(
      # draft names
      lambda_hat = lambda, p_hat = pc$p, b_hat = b, w_hat = pc$w, e_hat = pc$e,
      convergence = est$convergence,
      # canonical aliases
      coefficients = b, fitted.values = yhat, residuals = resid,
      e = pc$e, H_signal = H_signal, S = S,
      objective = est$value, value = est$value,
      converged = (est$convergence == 0L), method = "dual",
      # standard errors
      se_method = se_method, boot = boot, se_beta = NULL, vcov = NULL,
      # extras
      Z = Z, p0 = p0, v = v, w0 = w0, nu = nu, X = X, y = y, IV = IV,
      iv_scale = iv_scale, moment_resid = moment_resid,
      N = N, K = K, P = P, M = M, call = match.call()
    ),
    class = c("linreg_iv", "infometrics")
  )

  ## ---- standard errors (Z-estimator sandwich; Golan 2008, Sec. 3.3) --------
  if (se_method != "none") {
    V <- if (identical(se_method, "bootstrap"))
           .linreg_iv_boot_vcov(out, boot)
         else
           .linreg_iv_analytic_vcov(X, IVs, Z, v, nu, pc$p, pc$w, b, pc$e,
                                    resid - pc$e, meat = se_method)
    out$vcov <- V
    out$se_beta <- sqrt(diag(V))
  }
  out
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
vcov.linreg_iv <- function(object, type = c("sandwich", "delta"), ...) {
  if (missing(type) && !is.null(object$vcov)) return(object$vcov)
  type <- match.arg(type)
  IVs <- sweep(object$IV, 2L, object$iv_scale, "/")
  .linreg_iv_analytic_vcov(object$X, IVs, object$Z, object$v, object$nu,
                           object$p_hat, object$w_hat, object$coefficients,
                           object$e_hat, object$residuals - object$e_hat,
                           meat = type)
}

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
summary.linreg_iv <- function(object, se_method = NULL,
                              digits = max(3L, getOption("digits") - 3L), ...) {
  ident <- if (object$P == object$K) "just-identified" else "over-identified"
  meth <- if (is.null(se_method)) object$se_method
          else match.arg(se_method, c("sandwich", "delta", "bootstrap"))
  if (is.null(meth) || identical(meth, "none")) meth <- "sandwich"
  V <- if (identical(meth, object$se_method) && !is.null(object$vcov)) object$vcov
       else if (identical(meth, "bootstrap")) .linreg_iv_boot_vcov(object, object$boot)
       else vcov(object, type = meth)
  se <- sqrt(diag(V))
  z  <- object$coefficients / se
  pv <- 2 * stats::pnorm(-abs(z))
  tab <- cbind(Estimate = object$coefficients, `Std. Error` = se,
               `z value` = z, `Pr(>|z|)` = pv)

  cat("\nStochastic-moments GME-IV (Golan 2008, pp. 89-91)\n")
  cat("Call: "); print(object$call); cat("\n")
  cat(sprintf("  N=%d  K=%d  instruments=%d (%s)  M=%d  nu=%.3g\n",
              object$N, object$K, object$P, ident, object$M, object$nu))
  cat("\nResiduals (y - X beta):\n")
  print(summary(object$residuals, digits = digits))
  cat(sprintf("\nCoefficients (beta = Z p;  SE method: %s):\n", meth))
  stats::printCoefmat(tab, digits = digits, has.Pvalue = TRUE)
  cat(sprintf("\nnormalized signal entropy S = %.4f   max|moment resid| = %s\n",
              object$S, format(object$moment_resid, digits = digits)))
  cat(sprintf("convergence: %s\n", if (object$converged) "yes (0)" else "no"))
  cat("\nNote: asymptotic standard errors. The sampling distribution of beta is\n",
      "support-bounded and can be skewed, so Wald intervals are approximate with\n",
      "weak instruments or small n.\n", sep = "")
  invisible(object)
}
