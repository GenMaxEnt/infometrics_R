# inverse_noise.R
# Noisy (stochastic-moment) inverse problem via Generalized Maximum Entropy /
# Generalized Cross-Entropy, solved through the concentrated (dual) model.
# Formula-interface sibling of inverse_ce(): the response is the moment
# vector and each right-hand-side term is a state / support point, but the
# moments are treated as noisy (y = X p + e).
#
# References:
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145.
#     Section 6.1 (GME/GCE dual concentrated model).

# ---- internal: dual solver engine (reused by the bootstrap) -----------------

#' @keywords internal
.inverse_noise_engine <- function(X, y, v, nu, p0, w0, con,
                                  init = rep(0, nrow(X))) {
  Tn <- nrow(X)
  logp0 <- log(p0); logw0 <- log(w0)

  .signal <- function(lambda) {
    z <- logp0 + as.vector(lambda %*% X) / (1 - nu)
    zmax <- max(z); ez <- exp(z - zmax); s <- sum(ez)
    list(logOmega = zmax + log(s), p = ez / s)
  }
  .noise <- function(lambda) {
    a <- logw0 + outer(lambda, v) / nu
    amax <- apply(a, 1L, max)
    ea <- exp(a - amax); rs <- rowSums(ea)
    w <- ea / rs
    list(logPsi = amax + log(rs), w = w, e = as.vector(w %*% v))
  }
  obj <- function(lambda) {
    sg <- .signal(lambda); ns <- .noise(lambda)
    -sum(lambda * y) + (1 - nu) * sg$logOmega + nu * sum(ns$logPsi)
  }
  grad <- function(lambda) {
    sg <- .signal(lambda); ns <- .noise(lambda)
    as.vector(X %*% sg$p + ns$e - y)              # Xp + e - y
  }

  fit <- stats::optim(par = init, fn = obj, gr = grad, method = "BFGS",
                      control = list(maxit = con$maxit, reltol = con$reltol))
  lambda <- fit$par
  sg <- .signal(lambda); ns <- .noise(lambda)
  p <- sg$p; w <- ns$w; e <- ns$e
  Ep <- as.vector(X %*% p)

  ## analytic dual Hessian (positive definite thanks to the noise term)
  cov_sig   <- X %*% (p * t(X)) - tcrossprod(Ep)
  var_noise <- as.vector(w %*% (v^2)) - e^2
  hessian   <- (1 / (1 - nu)) * cov_sig + (1 / nu) * diag(var_noise, Tn)

  list(lambda = lambda, p = p, w = w, e = e, fitted = Ep, hessian = hessian,
       objective = fit$value, convergence = fit$convergence)
}

# ---- internal: standard errors (delta / sandwich / bootstrap) ---------------

#' @keywords internal
.inverse_noise_se <- function(X, y, p, resid, w, v, nu, p0, w0, hessian,
                              method, B, con) {
  Tn <- nrow(X); K <- ncol(X)
  Ep <- as.vector(X %*% p)
  ## Jacobian J_{k,t} = dp_k/dlambda_t = (1/(1-nu)) p_k (x_tk - E_p[x_t])
  J  <- (1 / (1 - nu)) * (p * (t(X) - matrix(Ep, K, Tn, byrow = TRUE)))   # K x T

  se_from_V <- function(V)
    list(vcov_lambda = V,
         se_lambda   = sqrt(pmax(diag(V), 0)),
         se_p        = sqrt(pmax(diag(J %*% V %*% t(J)), 0)))

  if (method == "delta") {
    Hi <- tryCatch(solve(hessian), error = function(e) NULL)
    if (is.null(Hi)) return(NULL)
    return(se_from_V(Hi))
  }
  if (method == "sandwich") {
    Hi <- tryCatch(solve(hessian), error = function(e) NULL)
    if (is.null(Hi)) return(NULL)
    ## robust (HC0) sandwich: bread H^-1, meat = diag(residuals^2) (fitted noise)
    return(se_from_V(Hi %*% diag(resid^2, Tn) %*% Hi))
  }
  if (method == "bootstrap") {
    ## residual bootstrap: rows fixed, resample the fitted noise, refit.
    ## Rows fixed => lambda* aligns to the original moments.
    Lb <- matrix(NA_real_, B, Tn); Pb <- matrix(NA_real_, B, K)
    for (b in seq_len(B)) {
      ys <- Ep + sample(resid, Tn, replace = TRUE)
      fb <- tryCatch(.inverse_noise_engine(X, ys, v, nu, p0, w0, con),
                     error = function(e) NULL)
      if (!is.null(fb)) { Lb[b, ] <- fb$lambda; Pb[b, ] <- fb$p }
    }
    Vlam <- stats::cov(Lb, use = "pairwise.complete.obs")
    return(list(vcov_lambda = Vlam,
                se_lambda   = apply(Lb, 2L, stats::sd, na.rm = TRUE),
                se_p        = apply(Pb, 2L, stats::sd, na.rm = TRUE)))
  }
  NULL
}

# ---- estimator --------------------------------------------------------------

#' Noisy Inverse Problem via Generalized Cross-Entropy (Formula Interface)
#'
#' Estimates a K-dimensional signal distribution \code{p} together with
#' per-moment noise distributions \code{w} for the stochastic-moment problem
#' \eqn{y = X p + e}, \eqn{e_t = \sum_j v_j w_{tj}}, by solving the
#' unconstrained dual (concentrated) model of Golan (2008, Section 6.1). This
#' is the noisy-moment counterpart of \code{\link{inverse_ce}}.
#'
#' Unlike \code{\link{inverse_ce}}, the moments are treated as \emph{noisy}:
#' the constraint is \eqn{y = X p + e} rather than \eqn{y = X p}. The estimator
#' therefore shrinks \code{p} toward its prior \code{p0}; how much depends on
#' the width of the noise support \code{v}. As \code{v} collapses toward 0 the
#' solution approaches the pure-CE solution of \code{\link{inverse_ce}}.
#'
#' @details
#' The dual objective is minimised here (convex) to match the package
#' convention used by \code{\link{inverse_ce}}:
#' \deqn{\ell(\lambda) = -\sum_t \lambda_t y_t + (1-\nu)\log\Omega(\lambda)
#'       + \nu \sum_t \log\Psi_t(\lambda),}
#' with signal partition function
#' \eqn{\Omega = \sum_k p_{0,k} \exp(\sum_t \lambda_t x_{tk}/(1-\nu))} and
#' per-moment noise partition functions
#' \eqn{\Psi_t = \sum_j w_{0,tj}\exp(\lambda_t v_j/\nu)}. The weight
#' \eqn{\nu \in (0,1)} splits entropy between signal and noise; \eqn{\nu = 0.5}
#' gives the standard equal-weight GME/GCE. Minimising this convex objective is
#' equivalent to maximising the joint signal-plus-noise entropy subject to the
#' noisy moment constraints. The estimates are recovered as
#' \eqn{\hat p_k \propto p_{0,k}\exp(\sum_t \hat\lambda_t x_{tk}/(1-\nu))} and
#' \eqn{\hat w_{tj} \propto w_{0,tj}\exp(\hat\lambda_t v_j/\nu)}, with
#' \eqn{\hat e_t = \sum_j v_j \hat w_{tj}}. Uniform priors (\code{p0 = NULL},
#' \code{w0 = NULL}) give the GME version; user-supplied priors give GCE.
#'
#' @section Standard errors:
#' Structurally this is a GME regression \eqn{y = Xp + e} with \code{p} on the
#' simplex, so the \eqn{T} moment rows are the effective sample size and sampling
#' SEs are meaningful. Because \eqn{\hat\lambda} solves
#' \eqn{Xp(\lambda)+e(\lambda)-y=0}, the implicit-function theorem gives
#' \eqn{\partial\hat\lambda/\partial y = H^{-1}}, so
#' \eqn{\mathrm{Var}(\hat\lambda) = H^{-1}\Sigma_e H^{-1}} and
#' \eqn{\mathrm{Var}(\hat p) = J\,\mathrm{Var}(\hat\lambda)\,J'} with
#' \eqn{J_{kt} = \partial p_k/\partial\lambda_t = (1-\nu)^{-1} p_k(x_{tk}-(Xp)_t)}.
#' Three methods are provided via \code{se_method} / \code{vcov(type=)} /
#' \code{summary(se_method=)}:
#' \describe{
#'   \item{\code{"sandwich"} (default)}{robust HC0
#'     \eqn{H^{-1}\,\mathrm{diag}(\hat e^2)\,H^{-1}} (meat = squared fitted noise).
#'     Accurate for \code{p} (validated against a Monte-Carlo SD and the
#'     bootstrap); per-element \code{lambda} SEs are noisier (each rests on one
#'     \eqn{\hat e_t^2}).}
#'   \item{\code{"delta"}}{the naive \eqn{H^{-1}}. Kept for comparison but it
#'     \strong{overstates the sampling SE by roughly ten-fold} (it is not a
#'     sampling covariance) -- do not use it for inference.}
#'   \item{\code{"bootstrap"}}{residual bootstrap: keep the rows fixed, resample
#'     the fitted noise \eqn{\hat e}, refit. Gives stable, aligned SEs for both
#'     \code{lambda} and \code{p}; the gold-standard cross-check.}
#' }
#' All three mildly understate (~15-20\%) because the finite-support entropy
#' penalty shrinks the fitted noise below the true noise dispersion (the
#' bootstrap shows the same, so it is a property of the estimator and of an
#' over-wide support, not a bug).
#'
#' @inheritSection inverse_ce Formula orientation
#'
#' @param formula A model formula. The response is the moment vector; each RHS
#'   term is a state. Usually written without an intercept, e.g.
#'   \code{y ~ x1 + x2 + x3 - 1}.
#' @param data A data frame (or environment) with one row per moment.
#' @param v Noise support. One of: \code{NULL} (default 3-point symmetric
#'   support of half-width \eqn{\max(3, \sqrt{2\log T})\,s} with
#'   \eqn{s = sd(y)}); a single whole number \code{M >= 2} (build an
#'   \code{M}-point symmetric support); or an explicit numeric vector of support
#'   points (should be symmetric around 0 for zero-mean noise). The default
#'   widens with \code{T} rather than using a fixed three-sigma rule, because
#'   the largest of \code{T} noise draws grows like \eqn{\sigma\sqrt{2\log T}};
#'   see the \emph{Support feasibility} section.
#' @param nu Entropy weight in \eqn{(0,1)} on the noise term; default 0.5
#'   (equal-weight GME/GCE).
#' @param p0 Prior signal probabilities, length K (= model-matrix columns);
#'   strictly positive, default uniform (GME). A non-uniform prior gives GCE.
#' @param w0 Prior noise probabilities, a \code{T x M} matrix (one row per
#'   moment); strictly positive rows summing to 1, default uniform. Its column
#'   count must match \code{length(v)}.
#' @param se_method Which standard errors to compute and store: \code{"sandwich"}
#'   (default, accurate), \code{"delta"} (naive \eqn{H^{-1}}; overstates ~10x),
#'   \code{"bootstrap"} (residual bootstrap; slower), or \code{"none"} (skip).
#'   See the \emph{Standard errors} section.
#' @param B Integer number of bootstrap resamples when
#'   \code{se_method = "bootstrap"} (or \code{vcov}/\code{summary} with
#'   \code{type}/\code{se_method = "bootstrap"}); default 1000.
#' @param subset Optional vector specifying a subset of moments (rows) to use,
#'   as in \code{\link[stats]{lm}}.
#' @param na.action A function indicating how to handle \code{NA}s in the model
#'   frame, as in \code{\link[stats]{lm}}.
#' @param control A named list merged over the defaults and passed to
#'   \code{\link[stats]{optim}} (BFGS): \code{maxit} (default 500) and
#'   \code{reltol} (default 1e-10).
#' @param ... Currently unused.
#'
#' @section Support feasibility:
#' The dual is bounded below only if \code{y} can be represented as
#' \eqn{Xp + e} with \eqn{p} on the simplex and every \eqn{e_t} inside the noise
#' support \code{v}. If it cannot, the dual is \strong{unbounded}: the
#' multipliers diverge, the softmax saturates, and \code{p_hat} collapses to a
#' vertex of the simplex (a 0/1 corner) even though \code{optim} reports
#' convergence. \code{inverse_noise()} therefore checks the first-order
#' condition \eqn{y - Xp - e = 0} at the optimum, reports it as
#' \code{foc_residual}, and \strong{warns} when it is not satisfied. The remedy
#' is a wider noise support \code{v}.
#'
#' @return An object of class \code{c("inverse_noise", "infometrics")}: a list
#'   with components
#'   \describe{
#'     \item{\code{p_hat}}{Named K-vector of signal probabilities.}
#'     \item{\code{foc_residual}}{\eqn{\max|y - Xp - e|} at the optimum; near 0
#'       for a healthy fit (see \emph{Support feasibility}).}
#'     \item{\code{lambda_hat}}{Named T-vector of Lagrange multipliers.}
#'     \item{\code{w_hat}}{\code{T x M} matrix of noise probabilities.}
#'     \item{\code{fitted.values}}{Signal-fitted moments \eqn{X p} (length T).}
#'     \item{\code{residuals}}{Moment residuals \eqn{y - X p}, equal to the
#'       estimated noise \eqn{e_t = \sum_j v_j w_{tj}} at the optimum.}
#'     \item{\code{hessian}}{Analytic T-by-T dual Hessian (positive definite).}
#'     \item{\code{vcov_lambda}, \code{se_lambda}, \code{se_p}}{Covariance of
#'       \eqn{\lambda} and standard errors of \eqn{\lambda}/\eqn{p} by
#'       \code{se_method} (\code{NULL} when \code{se_method = "none"}).}
#'     \item{\code{se_method}}{The SE method stored.}
#'     \item{\code{H_signal}, \code{H_error}, \code{H_error0}}{Shannon entropies
#'       of the signal \eqn{H(\hat p)}, the noise \eqn{H(\hat w)}, and the noise
#'       prior \eqn{H(w_0)}.}
#'     \item{\code{S}, \code{S_p}}{Normalized \emph{signal} entropy
#'       \eqn{H(\hat p)/H(p_0)} in [0, 1] (pseudo-R2 = 1 - S).}
#'     \item{\code{S_w}}{Normalized \emph{noise} entropy \eqn{H(\hat w)/H(w_0)}.}
#'     \item{\code{objective}}{Dual objective at the optimum.}
#'     \item{\code{nu}, \code{prior}, \code{noise_prior}, \code{support}}{The
#'       entropy weight, signal prior \code{p0}, noise prior \code{w0}, and the
#'       noise support \code{v}.}
#'     \item{\code{X}, \code{y}, \code{control}}{Stored for SE recomputation and
#'       the bootstrap.}
#'     \item{\code{convergence}}{Raw \code{optim} convergence code.}
#'     \item{\code{method}, \code{terms}, \code{model}, \code{call}}{Solver and
#'       standard model components.}
#'   }
#'
#' @references Golan, A. (2008). \emph{Information and Entropy Econometrics -
#'   A Review and Synthesis}. Foundations and Trends in Econometrics,
#'   2(1-2), 1-145. Section 6.1.
#'
#' @seealso \code{\link{inverse_ce}} for the exact-moment (pure) version;
#'   \code{\link{linreg}} for the GME/GCE regression estimator;
#'   \code{\link{fano_bounds}} for the Fano error bound on \code{p_hat}.
#'
#' @examples
#' set.seed(123)
#' X <- matrix(runif(30), nrow = 10, ncol = 3)   # T = 10 moments, K = 3 states
#' p_true <- c(0.5, 0.2, 0.3)
#' y <- as.vector(X %*% p_true) + rnorm(nrow(X), 0, 0.3)
#' dat <- data.frame(y = y, x1 = X[, 1], x2 = X[, 2], x3 = X[, 3])
#' fit <- inverse_noise(y ~ x1 + x2 + x3 - 1, data = dat)
#' fit
#' summary(fit)                 # p table (SE + t), signal/noise entropy, Fano
#' vcov(fit)                    # sandwich Cov(lambda)
#' fano_bounds(fit)             # Fano error bound for p_hat
#'
#' @importFrom stats model.frame model.matrix model.response optim sd cov printCoefmat vcov
#' @export
inverse_noise <- function(formula, data, v = NULL, nu = 0.5,
                          p0 = NULL, w0 = NULL,
                          se_method = c("sandwich", "delta", "bootstrap", "none"),
                          B = 1000L,
                          subset, na.action,
                          control = list(), ...) {

  se_method <- match.arg(se_method)

  ## ---- model frame (lm-style) --------------------------------------------
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- as.vector(stats::model.response(mf, "numeric"))
  X <- stats::model.matrix(mt, mf)
  if (attr(mt, "intercept") == 1L)
    warning("Formula includes an intercept; it becomes an extra 'state' with ",
            "constant moment 1. Use '- 1' to drop it if unintended.")

  Tn <- nrow(X)   # number of moments (T)
  K  <- ncol(X)   # number of states  (K)

  if (anyNA(y) || anyNA(X))
    stop("Missing values present after na.action.")
  if (length(y) != Tn)
    stop("length(y) (", length(y), ") must equal nrow(model matrix) (", Tn,
         ").")

  ## ---- nu and signal prior p0 --------------------------------------------
  if (length(nu) != 1L || !is.finite(nu) || nu <= 0 || nu >= 1)
    stop("'nu' must be a single value strictly between 0 and 1. Got nu = ", nu,
         ".")

  if (is.null(p0)) p0 <- rep(1 / K, K)
  p0 <- as.vector(p0)
  if (length(p0) != K)
    stop("length(p0) (", length(p0), ") must equal the number of states (", K,
         ").")
  if (any(p0 <= 0)) stop("All prior masses 'p0' must be strictly positive.")
  if (abs(sum(p0) - 1) > 1e-8) {
    warning("Prior 'p0' did not sum to 1; renormalizing."); p0 <- p0 / sum(p0)
  }

  ## ---- resolve noise support v and prior w0 ------------------------------
  sdY  <- stats::sd(y)
  span <- if (is.finite(sdY) && sdY > 0) 3 * sdY else 3
  ## The largest of Tn noise draws grows like sigma * sqrt(2 log Tn), so a
  ## fixed 3-sigma support becomes infeasible at large Tn -- y is then not
  ## representable as X p + e, the dual is unbounded below, lambda diverges
  ## and p collapses to a vertex of the simplex. Widen the span with Tn.
  espan <- max(3, sqrt(2 * log(max(Tn, 2L)))) / 3 * span
  is_count <- function(z) length(z) == 1L && is.finite(z) &&
                          abs(z - round(z)) < 1e-8

  if (is.null(v)) {
    M <- if (is.null(w0)) 3L else ncol(w0)
    v <- seq(-espan, espan, length.out = M)
  } else if (is_count(v)) {
    M <- as.integer(round(v))
    if (M < 2L)
      stop("A scalar 'v' is read as a support-point count and must be >= 2.")
    v <- seq(-espan, espan, length.out = M)
  } else {
    v <- as.vector(v)
    M <- length(v)
    if (M < 2L) stop("'v' must have at least 2 support points.")
  }
  if (is.null(w0)) {
    w0 <- matrix(1 / M, nrow = Tn, ncol = M)
  } else {
    w0 <- as.matrix(w0)
    if (ncol(w0) != M)
      stop("ncol(w0) (", ncol(w0), ") must equal length(v) (", M, ").")
    if (nrow(w0) != Tn)
      stop("nrow(w0) (", nrow(w0), ") must equal the number of moments (", Tn,
           ").")
    if (any(w0 <= 0))
      stop("All noise prior masses 'w0' must be strictly positive.")
    rs0 <- rowSums(w0)
    if (any(abs(rs0 - 1) > 1e-8)) {
      warning("Rows of 'w0' did not sum to 1; renormalizing."); w0 <- w0 / rs0
    }
  }
  if (abs(mean(v)) > 1e-8 * max(1, espan))
    warning("Noise support 'v' is not centered at 0; noise will not be ",
            "mean-zero.")

  ## ---- solve the dual (concentrated) model -------------------------------
  con <- list(maxit = 500L, reltol = 1e-10)
  con[names(control)] <- control
  eng <- .inverse_noise_engine(X, y, v, nu, p0, w0, con)
  if (eng$convergence != 0L)
    warning("optim did not converge (code ", eng$convergence, ").")

  lambda <- eng$lambda
  names(lambda) <- if (is.null(rownames(X)))
    paste0("moment", seq_len(Tn)) else rownames(X)
  p <- eng$p; names(p) <- colnames(X)
  w <- eng$w; e <- eng$e
  Ep <- eng$fitted; resid <- as.vector(y - Ep)
  hessian <- eng$hessian

  ## ---- feasibility / unbounded-dual check --------------------------------
  ## At the optimum the first-order condition is y - X p - e = 0. When the
  ## supports cannot represent y the dual is unbounded below: optim walks off
  ## to |lambda| -> Inf (still reporting convergence 0), the softmax saturates
  ## and p collapses to a vertex, leaving this FOC residual far from 0.
  foc <- max(abs(resid - e))
  if (!all(is.finite(lambda)) || max(abs(lambda)) > 1e8 ||
      foc > 1e-4 * max(1, sdY))
    warning("The GME/GCE dual appears unbounded: the noise support is too ",
            "narrow to represent y (max|y - X*p - e| = ",
            format(foc, digits = 3), "). Estimates are unreliable and 'p' may ",
            "collapse to a vertex of the simplex. Widen the noise support 'v'.")

  ## ---- normalized entropy: signal AND noise (mirror matrix_gce) ----------
  ent  <- function(mm) { z <- mm[mm > 0]; -sum(z * log(z)) }
  H_p  <- shannon_entropy(p); H_q <- shannon_entropy(p0)
  S    <- if (H_q > 0) H_p / H_q else NA_real_
  H_w  <- ent(w); H_w0 <- ent(w0)
  S_w  <- if (H_w0 > 0) H_w / H_w0 else NA_real_

  ## ---- standard errors ---------------------------------------------------
  se <- if (se_method == "none") NULL else
    .inverse_noise_se(X, y, p, resid, w, v, nu, p0, w0, hessian, se_method, B,
                      con)
  if (is.null(se)) se <- list(vcov_lambda = NULL, se_lambda = NULL, se_p = NULL)
  if (!is.null(se$se_lambda)) names(se$se_lambda) <- names(lambda)
  if (!is.null(se$se_p))      names(se$se_p)      <- names(p)
  if (!is.null(se$vcov_lambda))
    dimnames(se$vcov_lambda) <- list(names(lambda), names(lambda))

  structure(
    list(
      p_hat         = p,
      lambda_hat    = lambda,
      w_hat         = w,
      fitted.values = Ep,
      residuals     = resid,                   # = estimated noise e at optimum
      hessian       = hessian,
      vcov_lambda   = se$vcov_lambda,
      se_lambda     = se$se_lambda,
      se_p          = se$se_p,
      se_method     = se_method,
      H_signal      = H_p,                     # CLAUDE.md canonical field
      H_error       = H_w,
      H_error0      = H_w0,
      S             = S,                       # CLAUDE.md canonical (signal)
      S_p           = S,
      S_w           = S_w,
      objective     = eng$objective,           # CLAUDE.md canonical field
      foc_residual  = foc,                     # max|y - X p - e| at the optimum
      nu            = nu,
      prior         = p0,
      noise_prior   = w0,
      support       = v,
      X             = X,
      y             = y,
      control       = con,
      convergence   = eng$convergence,
      method        = "dual",                  # CLAUDE.md canonical field
      terms         = mt,
      model         = mf,
      call          = cl
    ),
    class = c("inverse_noise", "infometrics")
  )
}

## ----------------------------------------------------------------------------
## S3 methods (lm()-style)
## ----------------------------------------------------------------------------

#' @describeIn inverse_noise Estimated signal probabilities \code{p}.
#' @param object,x An \code{inverse_noise} object.
#' @export
coef.inverse_noise <- function(object, ...) object$p_hat

#' @describeIn inverse_noise Fitted (signal) moments \eqn{X p}.
#' @export
fitted.inverse_noise <- function(object, ...) object$fitted.values

#' @describeIn inverse_noise Residuals \eqn{y - X p} (the estimated noise).
#' @export
residuals.inverse_noise <- function(object, ...) object$residuals

#' @describeIn inverse_noise Covariance of the Lagrange multipliers by method.
#'   \code{type = "sandwich"} (default) is the robust
#'   \eqn{H^{-1}\mathrm{diag}(\hat e^2)H^{-1}}; \code{"delta"} is the naive
#'   \eqn{H^{-1}} (overstates ~10x); \code{"bootstrap"} refits \code{B} residual
#'   resamples. \strong{Note}: the default changed from the raw \eqn{H^{-1}} to
#'   the sandwich, since \eqn{H^{-1}} is not a valid sampling covariance here.
#' @param type Covariance method for \code{vcov}: \code{"sandwich"} (default),
#'   \code{"delta"}, or \code{"bootstrap"}.
#' @export
vcov.inverse_noise <- function(object, type = c("sandwich", "delta", "bootstrap"),
                               B = 1000L, ...) {
  type <- match.arg(type)
  if (!is.null(object$vcov_lambda) && identical(type, object$se_method)) {
    V <- object$vcov_lambda
  } else {
    se <- .inverse_noise_se(object$X, object$y, object$p_hat, object$residuals,
                            object$w_hat, object$support, object$nu,
                            object$prior, object$noise_prior, object$hessian,
                            type, B, object$control)
    V  <- if (is.null(se)) NULL else se$vcov_lambda
  }
  if (!is.null(V))
    dimnames(V) <- list(names(object$lambda_hat), names(object$lambda_hat))
  V
}

#' @describeIn inverse_noise Compact printed display.
#' @param digits Number of significant digits to print.
#' @export
print.inverse_noise <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Estimated probabilities (p):\n")
  print.default(format(x$p_hat, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\nLagrange multipliers (lambda):\n")
  print.default(format(x$lambda_hat, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\nnu = ", format(x$nu, digits = digits),
      ";  normalized entropy: signal S(p) = ", format(x$S, digits = digits),
      ", noise S(w) = ", format(x$S_w, digits = digits),
      "  (pseudo-R2 = ", format(1 - x$S, digits = digits), ")\n", sep = "")
  invisible(x)
}

#' @describeIn inverse_noise lm()-style summary; prints a \code{p} coefficient
#'   table with standard errors and t-stats, the signal/noise normalized
#'   entropies, and a Fano line. \code{se_method} selects the SE method
#'   (default: the one stored on the object).
#' @param se_method For \code{summary}: SE method to report (\code{"sandwich"},
#'   \code{"delta"}, or \code{"bootstrap"}); defaults to the stored method.
#' @export
summary.inverse_noise <- function(object,
                                  se_method = NULL, B = 1000L, ...) {
  meth <- if (is.null(se_method)) {
    if (is.null(object$se_method) || object$se_method == "none") "sandwich"
    else object$se_method
  } else match.arg(se_method, c("sandwich", "delta", "bootstrap"))

  if (!is.null(object$se_p) && identical(meth, object$se_method)) {
    se_p <- object$se_p; se_l <- object$se_lambda
  } else {
    se   <- .inverse_noise_se(object$X, object$y, object$p_hat, object$residuals,
                              object$w_hat, object$support, object$nu,
                              object$prior, object$noise_prior, object$hessian,
                              meth, B, object$control)
    se_p <- if (is.null(se)) NULL else se$se_p
    se_l <- if (is.null(se)) NULL else se$se_lambda
  }
  K  <- length(object$p_hat)
  fb <- attr(.fano_row_bounds(matrix(object$p_hat, 1L, K), K), "overall")

  structure(
    list(call = object$call, p = object$p_hat, se_p = se_p,
         lambda = object$lambda_hat, se_lambda = se_l,
         residuals = object$residuals, nu = object$nu,
         S = object$S, S_w = object$S_w, fano = fb, se_method = meth,
         convergence = object$convergence),
    class = "summary.inverse_noise"
  )
}

#' @describeIn inverse_noise Print method for the summary object.
#' @export
print.summary.inverse_noise <- function(x,
                                         digits = max(3L, getOption("digits") - 3L),
                                         ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Residuals / estimated noise (y - Xp = e):\n")
  if (length(x$residuals) > 5L) print(summary(x$residuals, digits = digits))
  else print.default(format(x$residuals, digits = digits), print.gap = 2L,
                     quote = FALSE)

  cat("\nCoefficients (estimated probabilities p; SE method: ", x$se_method,
      "):\n", sep = "")
  if (!is.null(x$se_p)) {
    ptab <- cbind(Estimate = x$p, `Std. Error` = x$se_p,
                  `t value` = x$p / x$se_p)
    stats::printCoefmat(ptab, digits = digits, has.Pvalue = FALSE)
  } else {
    print.default(cbind(Estimate = format(x$p, digits = digits)),
                  quote = FALSE, print.gap = 2L)
  }

  cat("\nLagrange multipliers (lambda):\n")
  if (!is.null(x$se_lambda)) {
    ltab <- cbind(Estimate = x$lambda, `Std. Error` = x$se_lambda,
                  `z value` = x$lambda / x$se_lambda)
    stats::printCoefmat(ltab, digits = digits, has.Pvalue = FALSE)
  } else {
    print.default(format(x$lambda, digits = digits), print.gap = 2L,
                  quote = FALSE)
  }

  cat("\nnu = ", format(x$nu, digits = digits),
      ";  normalized entropy: signal S(p) = ", format(x$S, digits = digits),
      ", noise S(w) = ", format(x$S_w, digits = digits),
      "   pseudo-R2 = ", format(1 - x$S, digits = digits), "\n", sep = "")
  cat("Fano (sec 7.5): modal error pe = ",
      format(x$fano["mean_pe"], digits = digits), " >= bound ",
      format(x$fano["mean_pe_lower"], digits = digits), "\n", sep = "")
  cat("Convergence code:", x$convergence,
      if (x$convergence == 0L) "(converged)" else "(see ?optim)", "\n")

  cat("\nNote: SE method '", x$se_method,
      "'. The sandwich (HC0) is the accurate default (validated vs a\n",
      "Monte-Carlo SD and the residual bootstrap); 'delta' (H^-1) overstates the\n",
      "sampling SE ~10x; per-element lambda SEs are noisiest under the sandwich\n",
      "(use se_method='bootstrap' for lambda). All SEs mildly understate (~15-20%)\n",
      "due to the finite-support shrinkage of the fitted noise.\n", sep = "")
  invisible(x)
}

## ----------------------------------------------------------------------------
## Fano error bounds for the recovered distribution p  (Golan 2008, sec 7.5)
## ----------------------------------------------------------------------------

#' Fano Error Bounds for a Noisy Inverse-Problem (GME/GCE) Fit
#'
#' The recovered \code{p_hat} is a single distribution over the \eqn{K} states, so
#' Golan's (2008, sec 7.5) Fano bound applies to it directly (a parallel of
#' \code{\link{fano_bounds.inverse_ce}}): a modal classifier that guesses
#' \eqn{\arg\max_k p_k} errs with probability \eqn{pe = 1 - \max_k p_k}, and
#' \eqn{pe \ge S(p) - \log 2/\log K} with \eqn{S(p) = H(p)/\log K}
#' (uniform-reference, distinct from the prior-relative \code{S}/\code{S_w} the
#' object reports). An information-theoretic error bound, \emph{not} a sampling SE
#' (for those see \code{\link{vcov.inverse_noise}} / \code{summary}).
#'
#' @param object An \code{inverse_noise} object.
#' @param ... Unused.
#' @return A one-row data frame with columns \code{p_max}, \code{pe}, \code{H},
#'   \code{S}, and \code{pe_lower}; an \code{"overall"} attribute mirrors the row.
#' @references Golan, A. (2008). Information and Entropy Econometrics.
#'   \emph{Foundations and Trends in Econometrics}, \strong{2}(1-2), 1-145
#'   (sec 7.5); Fano, R. (1961). \emph{Transmission of Information}.
#' @seealso \code{\link{fano_bounds}}, \code{\link{inverse_noise}}
#' @export
fano_bounds.inverse_noise <- function(object, ...) {
  p <- object$p_hat
  .fano_row_bounds(matrix(p, 1L, length(p)), length(p))
}
