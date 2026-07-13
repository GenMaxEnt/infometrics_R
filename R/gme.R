# ============================================================
# gme.R — Generalized Maximum Entropy / Cross-Entropy estimator
#         for linear regression (Phase 2)
# ============================================================

# ------------------------------------------------------------------
# Internal: dual concentrated solver
# ------------------------------------------------------------------

#' @keywords internal
.gme_dual <- function(y, X, Z, p0, V, w0, nu, con) {
  n <- length(y)

  # log partition function for each coefficient k (K-vector)
  # Omega_k(lambda) = sum_m p0_km exp(-z_km * score_k / (1-nu))
  # where score_k = sum_t lambda_t x_tk = (lambda %*% X)[k]
  .log_omega_k <- function(lambda) {
    scores <- drop(lambda %*% X)                      # K-vector
    e_km   <- -sweep(Z, 1, scores, "*") / (1 - nu)   # K x M
    m_k    <- apply(e_km, 1, max)                     # row maxima (log-sum-exp shift)
    log(rowSums(p0 * exp(e_km - m_k))) + m_k
  }

  # log partition function for each observation t (n-vector)
  # Psi_t(lambda) = sum_j w0_tj exp(-lambda_t v_j / nu)
  .log_psi_t <- function(lambda) {
    e_tj <- -outer(lambda, V) / nu   # n x J
    m_t  <- apply(e_tj, 1, max)     # row maxima (log-sum-exp shift)
    log(rowSums(w0 * exp(e_tj - m_t))) + m_t
  }

  # dual objective to minimise (Golan 2008, below Eq. 6.9, nu-weighted form)
  # min_lambda { sum(lambda*y) + (1-nu)*sum(log Omega_k) + nu*sum(log Psi_t) }
  .gme_obj <- function(lambda) {
    sum(lambda * y) +
      (1 - nu) * sum(.log_omega_k(lambda)) +
      nu       * sum(.log_psi_t(lambda))
  }

  # analytic gradient: d obj / d lambda_s = y_s - (X beta)_s - e_s
  # (the (1-nu) and nu scalings cancel in the derivative)
  .gme_grad <- function(lambda) {
    scores <- drop(lambda %*% X)
    e_km   <- -sweep(Z, 1, scores, "*") / (1 - nu)
    m_k    <- apply(e_km, 1, max)
    p_km   <- p0 * exp(e_km - m_k)
    p_km   <- p_km / rowSums(p_km)
    beta   <- rowSums(Z * p_km)          # K-vector: E_p[Z]

    e_tj   <- -outer(lambda, V) / nu
    m_t    <- apply(e_tj, 1, max)
    w_tj   <- w0 * exp(e_tj - m_t)
    w_tj   <- w_tj / rowSums(w_tj)
    e_t    <- drop(w_tj %*% V)           # n-vector: E_w[V]

    as.vector(y - X %*% beta - e_t)
  }

  stats::optim(
    par     = rep(0, n),
    fn      = .gme_obj,
    gr      = .gme_grad,
    method  = "BFGS",
    control = list(
      maxit  = con$maxit,
      reltol = con$tol,
      trace  = as.integer(con$trace)
    )
  )
}


# ------------------------------------------------------------------
# Exported: gme()
# ------------------------------------------------------------------

#' Generalized Maximum Entropy / Cross-Entropy Regression
#'
#' @description
#' Estimates a linear regression model using the Generalized Maximum Entropy
#' (GME) or Generalized Cross-Entropy (GCE) estimator of Golan, Judge and
#' Miller (1996), following the unified treatment in Golan (2008, Chapter 6).
#'
#' The model is \eqn{y = X\hat{\beta} + \hat{\varepsilon}}, where each
#' coefficient \eqn{\hat{\beta}_k = \sum_m z_{km} \hat{p}_{km}} is the
#' expected value of a probability distribution over a signal support \eqn{Z_k},
#' and each residual \eqn{\hat{\varepsilon}_t = \sum_j v_j \hat{w}_{tj}} is the
#' expected value over an error support \eqn{V}.
#'
#' When \code{p0} and \code{w0} are \code{NULL} (uniform priors), the estimator
#' maximises \eqn{H(p) + H(w)} and reduces to GME. With user-supplied priors it
#' minimises \eqn{D(p \| p_0) + D(w \| w_0)} and becomes GCE. The relationship
#' mirrors \code{me()}: ME is CE with uniform prior; GME is GCE with uniform
#' priors.
#'
#' @details
#' The dual concentrated problem minimised by BFGS (Golan 2008, below Eq. 6.9,
#' nu-weighted form) is:
#' \deqn{
#'   \min_\lambda \Bigl\{
#'     \sum_t \lambda_t y_t
#'     + (1-\nu)\sum_k \log\Omega_k(\lambda)
#'     + \nu\sum_t \log\Psi_t(\lambda)
#'   \Bigr\}
#' }
#' where
#' \deqn{\Omega_k(\lambda) = \sum_m p_{0,km}
#'   \exp\!\Bigl(-z_{km}\tfrac{\sum_t \lambda_t x_{tk}}{1-\nu}\Bigr),}
#' \deqn{\Psi_t(\lambda) = \sum_j w_{0,tj}
#'   \exp\!\Bigl(-\tfrac{\lambda_t v_j}{\nu}\Bigr).}
#' Setting \eqn{\nu = 0.5} recovers standard GME (Golan 2008, Eq. 6.5) because
#' maximising \eqn{0.5\,H(p) + 0.5\,H(w)} has the same argmax as
#' \eqn{H(p) + H(w)}.
#'
#' Both partition functions are computed with the log-sum-exp trick for
#' numerical stability (CLAUDE.md).
#'
#' @param formula An object of class \code{\link{formula}}: a symbolic
#'   description of the model to be fitted, e.g. \code{y ~ x1 + x2}.
#' @param data A data frame containing the variables in \code{formula}.
#' @param Z K x M numeric matrix of signal support points, one row per
#'   coefficient. Default: data-driven via \code{\link{default_supports}}.
#' @param p0 K x M numeric matrix of signal prior probabilities (rows must sum
#'   to 1). Default: uniform (\code{1/M}).  Supplying \code{p0} activates GCE.
#' @param V Numeric vector of length J giving the error support (must be
#'   symmetric around zero). Default: data-driven via
#'   \code{\link{default_supports}}.
#' @param w0 n x J numeric matrix of error prior probabilities (rows must sum
#'   to 1). Default: uniform (\code{1/J}).  Supplying \code{w0} activates GCE.
#' @param nu Single numeric in \eqn{(0, 1)}: weight assigned to error entropy
#'   relative to signal entropy. \code{nu = 0.5} (default) gives standard GME.
#' @param method Character: \code{"dual"} (default) solves the concentrated
#'   dual via BFGS. \code{"primal"} is not yet implemented.
#' @param control Named list of solver options: \code{maxit} (default 500L),
#'   \code{tol} (default 1e-10), \code{trace} (default \code{FALSE}).
#'
#' @return An object of class \code{c("infometrics_gme", "infometrics")} with
#'   components:
#'   \describe{
#'     \item{\code{coefficients}}{Named K-vector of estimated regression
#'       coefficients \eqn{\hat{\beta}_k = \sum_m z_{km} \hat{p}_{km}}.}
#'     \item{\code{residuals}}{n-vector of estimated errors
#'       \eqn{\hat{\varepsilon}_t = \sum_j v_j \hat{w}_{tj}}.}
#'     \item{\code{fitted.values}}{n-vector \eqn{\hat{y} = X\hat{\beta}}.}
#'     \item{\code{vcov}}{K x K estimated variance-covariance matrix
#'       (Golan 2008, p. 96).}
#'     \item{\code{lambda_hat}}{n-vector of optimal Lagrange multipliers.}
#'     \item{\code{p_hat}}{K x M matrix of estimated signal probabilities.}
#'     \item{\code{w_hat}}{n x J matrix of estimated error probabilities.}
#'     \item{\code{H_signal}}{K-vector of per-coefficient Shannon entropy.}
#'     \item{\code{H_error}}{n-vector of per-observation error entropy.}
#'     \item{\code{S_p}}{K-vector of normalised signal entropy (per coeff).}
#'     \item{\code{S_P}}{Scalar overall normalised entropy.}
#'     \item{\code{pseudo_R2}}{Scalar \eqn{1 - S_P} (Golan 2008, Section 7.5).}
#'     \item{\code{sigma}}{Estimated residual scale.}
#'     \item{\code{converged}}{Logical: did the BFGS solver converge?}
#'     \item{\code{rank}, \code{df.residual}}{Integers n and n - K.}
#'   }
#'
#' @references
#' Golan, A. (2008). Information and Entropy Econometrics — A Review and
#' Synthesis. \emph{Foundations and Trends in Econometrics}, \strong{2}(1-2),
#' 1-145.
#'
#' Golan, A., Judge, G. and Miller, D. (1996). \emph{Maximum Entropy
#' Econometrics}. Wiley.
#'
#' @seealso \code{\link{me}} for the ME/CE estimator without regression
#'   structure; \code{\link{default_supports}} for automatic support
#'   construction.
#'
#' @examples
#' set.seed(42)
#' n  <- 40
#' x  <- rnorm(n)
#' y  <- 1 + 2 * x + rnorm(n, sd = 0.5)
#' df <- data.frame(y = y, x = x)
#'
#' fit <- gme(y ~ x, data = df)
#' print(fit)
#' coef(fit)
#' summary(fit)
#'
#' @export
gme <- function(formula, data, Z = NULL, p0 = NULL, V = NULL, w0 = NULL,
                nu = 0.5, method = c("dual", "primal"), control = list()) {

  mc     <- match.call()
  mf     <- model.frame(formula, data)
  mt     <- attr(mf, "terms")
  y      <- model.response(mf)
  X      <- model.matrix(mt, mf)

  n <- nrow(X)
  K <- ncol(X)

  if (!is.numeric(y) || !is.numeric(X))
    stop("y and X must be numeric.")

  # validate nu
  if (!is.numeric(nu) || length(nu) != 1L || nu <= 0 || nu >= 1)
    stop("nu must be a single number in (0, 1). Got nu = ", nu, ".")

  method <- match.arg(method)
  if (method == "primal")
    stop("method = \"primal\" is not yet implemented for gme().")

  # control defaults
  con <- list(maxit = 500L, tol = 1e-10, trace = FALSE)
  con[names(control)] <- control

  # determine support sizes early (needed for default_supports call)
  M_sig <- if (!is.null(p0)) ncol(p0) else if (!is.null(Z)) ncol(Z) else 5L
  J_err <- if (!is.null(w0)) ncol(w0) else if (!is.null(V)) length(V) else 3L

  # data-driven defaults (OLS-based) via default_supports()
  if (is.null(Z) || is.null(V)) {
    sd_y <- stats::sd(y)
    sp   <- default_supports(y, X, M_signal = M_sig, M_error = J_err)
    if (is.null(Z)) Z <- sp$Z
    if (is.null(V)) V <- sp$V
  }

  M <- ncol(Z)
  J <- length(V)

  if (is.null(p0)) p0 <- matrix(1 / M, nrow = K, ncol = M)
  if (is.null(w0)) w0 <- matrix(1 / J, nrow = n, ncol = J)

  # ---- input validation ------------------------------------------------
  if (nrow(X) != length(y))
    stop("nrow(X) must equal length(y). Got nrow(X) = ", nrow(X),
         ", length(y) = ", length(y), ".")
  if (!is.matrix(Z) || nrow(Z) != K || ncol(Z) != M)
    stop("Z must be a ", K, " x ", M, " matrix. Got dim(Z) = (",
         paste(dim(Z), collapse = ", "), ").")
  if (!is.matrix(p0) || nrow(p0) != K || ncol(p0) != M)
    stop("p0 must be a ", K, " x ", M, " matrix. Got dim(p0) = (",
         paste(dim(p0), collapse = ", "), ").")
  if (!is.numeric(V) || length(V) != J || length(V) < 2L)
    stop("V must be a numeric vector of length >= 2. Got length(V) = ",
         length(V), ".")
  if (!is.matrix(w0) || nrow(w0) != n || ncol(w0) != J)
    stop("w0 must be a ", n, " x ", J, " matrix. Got dim(w0) = (",
         paste(dim(w0), collapse = ", "), ").")
  if (abs(sum(V)) > 1e-8 * J)
    stop("V must be symmetric around zero (sum(V) must be ~0). Got sum(V) = ",
         sum(V), ".")
  if (any(p0 <= 0) || any(abs(rowSums(p0) - 1) > 1e-8))
    stop("p0 must have strictly positive entries with rows summing to 1.")
  if (any(w0 <= 0) || any(abs(rowSums(w0) - 1) > 1e-8))
    stop("w0 must have strictly positive entries with rows summing to 1.")

  # ---- solve dual -------------------------------------------------------
  opt    <- .gme_dual(y, X, Z, p0, V, w0, nu, con)
  lambda <- opt$par

  # ---- recover probabilities (log-sum-exp stabilised) -------------------
  scores <- drop(lambda %*% X)                     # K-vector
  e_km   <- -sweep(Z, 1, scores, "*") / (1 - nu)  # K x M
  m_k    <- apply(e_km, 1, max)
  p_km   <- p0 * exp(e_km - m_k)
  p_km   <- p_km / rowSums(p_km)

  e_tj   <- -outer(lambda, V) / nu                 # n x J
  m_t    <- apply(e_tj, 1, max)
  w_tm   <- w0 * exp(e_tj - m_t)
  w_tm   <- w_tm / rowSums(w_tm)

  # ---- derived quantities -----------------------------------------------
  beta  <- rowSums(Z * p_km)           # K-vector: beta_k = E_p[Z_k]
  e_hat <- drop(w_tm %*% V)            # n-vector: eps_t = E_w[V_t]
  y_hat <- drop(X %*% beta)            # n-vector: y_hat = X beta
  names(beta) <- colnames(X)

  H_U     <- rowSums(-p_km * log(p_km))   # K-vector: per-coeff signal entropy
  H_err   <- rowSums(-w_tm * log(w_tm))   # n-vector: per-obs error entropy
  H_R     <- rep(log(M), K)               # K-vector: uniform reference entropy
  S_p     <- H_U / log(M)                 # K-vector: normalised signal entropy
  S_P     <- mean(S_p)                    # scalar: overall normalised entropy
  pseudo_R2 <- 1 - S_P                    # Golan (2008) Section 7.5

  # variance (Golan 2008, p. 96)
  s2     <- sum(lambda^2) / n
  var_et <- drop(w_tm %*% V^2) - e_hat^2  # Var_t(V) = E[V^2] - E[V]^2
  omega2 <- mean(1 / var_et)^2
  vcov_m <- tryCatch(
    (s2 / omega2) * solve(crossprod(X)),
    error = function(e) {
      warning("vcov computation failed (possibly singular X'X). Returning NA matrix.")
      matrix(NA_real_, K, K)
    }
  )
  rownames(vcov_m) <- colnames(vcov_m) <- colnames(X)

  # ---- return object ----------------------------------------------------
  structure(
    list(
      call          = mc,
      terms         = mt,
      coefficients  = beta,
      residuals     = e_hat,
      fitted.values = y_hat,
      rank          = K,
      df.residual   = n - K,
      sigma         = sqrt(s2),
      vcov          = vcov_m,
      lambda_hat    = lambda,
      p_hat         = p_km,
      w_hat         = w_tm,
      H_signal      = H_U,
      H_error       = H_err,
      H_R           = H_R,
      S_p           = S_p,
      S_P           = S_P,
      pseudo_R2     = pseudo_R2,
      objective     = opt$value,
      converged     = (opt$convergence == 0L),
      method        = method,
      nu            = nu,
      Z             = Z,
      V             = V,
      model         = mf,
      xlevels       = .getXlevels(mt, mf),
      contrasts     = attr(X, "contrasts")
    ),
    class = c("infometrics_gme", "infometrics")
  )
}


# ------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------

#' @export
print.infometrics_gme <- function(x, digits = 4L, ...) {
  K <- x$rank
  M <- ncol(x$p_hat)
  J <- ncol(x$w_hat)
  n <- nrow(x$w_hat)

  cat("\nGeneralized Maximum Entropy Regression\n")
  cat("---------------------------------------\n")
  cat("n =", n, "  K =", K, "  M =", M, "  J =", J,
      "  nu =", x$nu, "\n")
  cat("Method:", x$method,
      " Converged:", x$converged, "\n\n")
  cat("Coefficients:\n")
  print(round(x$coefficients, digits))
  cat("\nNormalised entropy (S_P): ", round(x$S_P,       digits), "\n")
  cat("Pseudo-R2:                ", round(x$pseudo_R2,  digits), "\n")
  invisible(x)
}

#' @export
summary.infometrics_gme <- function(object, digits = 4L, ...) {
  se    <- sqrt(diag(object$vcov))
  tstat <- object$coefficients / se
  tab   <- cbind(
    Estimate = object$coefficients,
    Std.Err  = se,
    t.value  = tstat,
    S_p      = object$S_p
  )

  cat("\nGeneralized Maximum Entropy Regression\n")
  cat("Call: "); print(object$call); cat("\n")

  cat("Residuals:\n")
  print(round(stats::quantile(object$residuals), digits))
  cat("\n")

  cat("Coefficients:\n")
  print(round(tab, digits))
  cat("\n")

  cat("Residual scale (sigma):  ", round(object$sigma,     digits), "\n")
  cat("Overall S_P:             ", round(object$S_P,       digits), "\n")
  cat("Pseudo-R2:               ", round(object$pseudo_R2, digits), "\n")
  cat("nu =", object$nu,
      " Converged:", object$converged, "\n")
  invisible(object)
}

#' @export
coef.infometrics_gme <- function(object, ...) {
  object$coefficients
}

#' @export
fitted.infometrics_gme <- function(object, ...) {
  object$fitted.values
}

#' @export
residuals.infometrics_gme <- function(object, ...) {
  object$residuals
}

#' @export
vcov.infometrics_gme <- function(object, ...) {
  object$vcov
}
