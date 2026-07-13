# panel_gce.R
# Dual GCE estimation of the one-way error-components panel regression model
#   y_nt = x_nt' beta + mu_n + eps_nt
# following Lee & Cheon (2014, CSAM 21(5)), eq. (3.12), extended with priors
# (p0, g0, w0) and the signal/noise weight nu. Panel/individual-effects sibling
# of gme()/linreg(): three reparameterized blocks solved by the concentrated
# dual over the N*T Lagrange multipliers lambda.
#
#   beta_k = sum_m z_km p_km,   p_km  prop to p0_km  exp( z_km (X'lambda)_k / nu )
#   mu_n   = sum_r f_nr g_nr,   g_nr  prop to g0_nr  exp( f_nr (sum_t lambda_nt) / nu )
#   e_nt   = sum_j v_j w_ntj,   w_ntj prop to w0_ntj exp( v_j lambda_nt / (1-nu) )
#   objective (MAXIMISED, fnscale = -1):
#     y'lambda - nu sum_k logOmega_k - nu sum_n logPhi_n - (1-nu) sum_nt logPsi_nt
#   gradient / FOC: y - X beta - rep(mu) - e = 0  (the model equation)
#
# The maximise form is the lambda -> -lambda equivalent of the paper's
# minimisation (eq. 3.12 uses exp(-.)), kept by user choice (with markov_*,
# multinomial_gce, linreg_iv). Standard errors follow Lee & Cheon (2014) Sec.
# 3.3: beta via the asymptotic variance Q^{-1} Xi Q^{-1} (eq. 3.17), with p and
# mu obtained by the delta method (see the SE block below).
#
# References:
#   Lee, S. and Cheon, S. (2014). Dual generalized maximum entropy estimation
#     for panel data regression models. Communications for Statistical
#     Applications and Methods, 21(5), 395-409.
#   Golan, A., Judge, G. and Perloff, J.M. (1996). A maximum entropy approach to
#     recovering information from multinomial response data. JASA, 91, 841-853.

#' Dual GCE Estimator for the One-Way Error-Components Panel Model
#'
#' @description
#' Estimates the panel regression \eqn{y_{nt} = x_{nt}'\beta + \mu_n +
#' \varepsilon_{nt}} by dual Generalized Cross-Entropy (Lee & Cheon 2014, eq.
#' 3.12), extended with priors \code{p0}, \code{g0}, \code{w0} and a
#' signal/noise weight \code{nu}. The coefficients \eqn{\beta}, the individual
#' effects \eqn{\mu_n}, and the errors \eqn{\varepsilon_{nt}} are each
#' reparameterized on bounded supports and recovered jointly via the
#' concentrated dual over the \eqn{N T} Lagrange multipliers.
#'
#' @details
#' With \eqn{\beta_k = \sum_m z_{km} p_{km}} (support \code{Z}, prior \code{p0}),
#' \eqn{\mu_n = \sum_r f_{nr} g_{nr}} (support \code{FF}, prior \code{g0}), and
#' \eqn{e_{nt} = \sum_j v_j w_{ntj}} (support \code{v}, prior \code{w0}), the dual
#' is maximized (\code{fnscale = -1}); \eqn{\beta} and \eqn{\mu} are weighted by
#' \code{nu}, the errors by \code{1 - nu}. The gradient is the model equation
#' \eqn{y - X\beta - \mu - e}, so at the optimum it holds to numerical zero
#' (\code{foc_residual}). As the supports widen the estimate approaches the
#' within (fixed-effects) estimate.
#'
#' Standard errors follow Lee & Cheon (2014, Sec. 3.3). Because \eqn{\mu} is
#' estimated, \eqn{\hat\beta} is identified from within-unit variation, so its
#' sampling variance is the within (fixed-effects) form \eqn{\hat\sigma^2_e
#' (\tilde X'\tilde X)^{-1}} (option \code{vcov_type = "within"}, the default and
#' the most accurate in Monte Carlo). The paper's stated asymptotic variance
#' \eqn{Q^{-1}\Xi Q^{-1}} (eq. 3.17) with the composite disturbance is available
#' as \code{vcov_type = "disturbance"} (it is conservative here). SEs for
#' \eqn{p} and \eqn{\mu} are obtained by the delta method from \eqn{Var(\hat
#' \beta)}: \eqn{SE(\hat p_{km}) = |p_{km}(z_{km}-\beta_k)| / Var_p(z_k) \cdot
#' SE(\hat\beta_k)} and \eqn{SE(\hat\mu_n) = \sqrt{\hat\sigma^2_e / T_n +
#' \bar x_n' Var(\hat\beta)\, \bar x_n}} (the latter a finite-\eqn{T} prediction
#' SE, since individual effects are not \eqn{\sqrt N}-consistent). \code{X} must
#' not include an intercept column (the individual effects subsume the level;
#' a constant column makes the within design singular).
#'
#' @param y Numeric response vector of length \code{N*tt}.
#' @param X \code{(N*tt)}-by-K design matrix.
#' @param Z Coefficient (signal) support: a K-by-M matrix, row k the support for
#'   \eqn{\beta_k}.
#' @param tt Integer number of time periods (balanced panel; \code{nrow(X)} must
#'   be a multiple of \code{tt}). The number of units is \code{N = nrow(X)/tt}.
#' @param FF Individual-effects support: an N-by-R matrix, row n the support for
#'   \eqn{\mu_n} (typically the same symmetric grid for every unit).
#' @param p0 Optional K-by-M signal prior (rows > 0, summing to 1); default
#'   uniform.
#' @param g0 Optional N-by-R effects prior (rows > 0, summing to 1); default
#'   uniform.
#' @param v Optional error support: \code{NULL} (default 3-point symmetric grid
#'   on \eqn{\pm 3} within-residual sd), a single whole number giving the number
#'   of support points, or an explicit numeric vector.
#' @param w0 Optional \code{(N*tt)}-by-J error prior (rows > 0, summing to 1);
#'   default uniform.
#' @param nu Signal/noise entropy weight in \eqn{(0,1)}; default 0.5.
#' @param layout How \code{lambda}/rows map to \eqn{(n, t)}: \code{"period"} (all
#'   units for t=1, then t=2, ...) or \code{"unit"} (all periods for unit 1, then
#'   unit 2, ...). Must match the ordering of \code{y}/\code{X}.
#' @param vcov_type Covariance estimator for \eqn{\hat\beta} (Lee & Cheon 2014,
#'   Sec. 3.3); one of:
#'   \describe{
#'     \item{\code{"within"}}{(default) the within/fixed-effects form
#'       \eqn{\hat\sigma^2_e (\tilde X'\tilde X)^{-1}} with \eqn{\tilde X} the
#'       within-demeaned design and \eqn{\hat\sigma^2_e = \sum \hat e^2 /
#'       (NT - N - K)}. The most accurate SE when \eqn{\mu} is well identified
#'       (\eqn{\hat\beta} is within-identified); matched simulation sampling SDs
#'       to within ~1\%.}
#'     \item{\code{"cluster"}}{the within estimator with a unit-clustered meat
#'       \eqn{(\tilde X'\tilde X)^{-1} [\sum_n \tilde X_n'\hat e_n \hat e_n'
#'       \tilde X_n] (\tilde X'\tilde X)^{-1}}; robust to heteroskedasticity and
#'       serial correlation in \eqn{\varepsilon}.}
#'     \item{\code{"disturbance"}}{the paper's stated asymptotic variance
#'       \eqn{Q^{-1}\Xi Q^{-1}} (eq. 3.17), i.e.
#'       \eqn{(X'X)^{-1} [\sum_n X_n'\hat u_n \hat u_n' X_n] (X'X)^{-1}} using the
#'       composite disturbance \eqn{\hat u = y - X\hat\beta} (unit-clustered).
#'       Conservative relative to \code{"within"} because it counts the
#'       between-unit variation that \eqn{\hat\mu} absorbs.}
#'   }
#'   The \eqn{p} and \eqn{\mu} SEs are delta-method transforms of the chosen
#'   \eqn{Var(\hat\beta)} (see \strong{Details}).
#' @param control Named list merged over defaults and passed to
#'   \code{\link[stats]{optim}} (BFGS): \code{maxit} (default 1000) and
#'   \code{reltol} (default 1e-12). \code{fnscale} is forced to -1.
#'
#' @return An object of class \code{c("panel_gce", "infometrics")} with
#'   \code{coefficients}/\code{b_hat} (\eqn{\beta}, K; \code{coef()} returns it),
#'   \code{mu_hat} (N individual effects), \code{e_hat} (N*tt errors),
#'   \code{p_hat}/\code{g_hat}/\code{w_hat} (reparameterization weights),
#'   \code{lambda_hat}, \code{fitted.values}, \code{residuals},
#'   \code{foc_residual}, \code{H_p}/\code{H_g}/\code{H_w}, \code{H_signal},
#'   \code{S_p}/\code{S}, \code{objective}/\code{value}, \code{converged}/
#'   \code{convergence}, \code{method}, \code{nu}, \code{v}, \code{layout}, the
#'   stored inputs, and \code{call}. Standard errors: \code{vcov} (K-by-K,
#'   accessible via \code{\link[stats]{vcov}}), \code{se_beta} (K), \code{se_p}
#'   (K-by-M), \code{se_mu} (N), \code{sigma2_eps}, and \code{vcov_type}.
#'
#' @references
#' Lee, S. and Cheon, S. (2014). Dual generalized maximum entropy estimation for
#' panel data regression models. \emph{Communications for Statistical
#' Applications and Methods}, \strong{21}(5), 395-409.
#'
#' @seealso \code{\link{gme}}, \code{\link{linreg}} for the (non-panel) GME/GCE
#'   regression.
#'
#' @examples
#' set.seed(1)
#' N <- 20L; T <- 4L; K <- 2L
#' mu   <- rnorm(N)                          # individual effects
#' id   <- rep(seq_len(N), each = T)         # unit-major ordering
#' X    <- matrix(rnorm(N * T * K), N * T, K)
#' y    <- as.vector(X %*% c(1.5, -0.8)) + mu[id] + rnorm(N * T, 0, 0.5)
#' Z  <- matrix(c(-10, 0, 10), nrow = K, ncol = 3, byrow = TRUE)  # beta support
#' FF <- matrix(c(-6, 0, 6),  nrow = N, ncol = 3, byrow = TRUE)   # mu support
#' fit <- panel_gce(y, X, Z, tt = T, FF = FF, layout = "unit")
#' coef(fit)                                 # ~ (1.5, -0.8) (within/FE estimate)
#'
#' @importFrom stats optim sd lm.fit ave residuals pnorm printCoefmat
#' @export
panel_gce <- function(y, X, Z, tt, FF,
                      p0 = NULL, g0 = NULL, v = NULL, w0 = NULL, nu = 0.5,
                      layout = c("period", "unit"),
                      vcov_type = c("within", "cluster", "disturbance"),
                      control = list()) {

  layout    <- match.arg(layout)
  vcov_type <- match.arg(vcov_type)

  ## ---- dimensions and validation -------------------------------------------
  X <- as.matrix(X); y <- as.numeric(y)
  nt <- nrow(X); kk <- ncol(X)
  if (length(y) != nt) stop("length(y) must equal nrow(X).")
  if (anyNA(y) || anyNA(X)) stop("'y' and 'X' must not contain NA.")
  if (nt %% tt != 0L) stop("nrow(X) must be a multiple of tt (balanced panel).")
  nn <- nt %/% tt

  ## unit index per row (used for the default-v scale and the SE block)
  id <- if (layout == "period") rep(seq_len(nn), times = tt) else rep(seq_len(nn), each = tt)

  Z <- as.matrix(Z)
  if (nrow(Z) != kk) stop("nrow(Z) must equal ncol(X) (one support row per beta_k).")
  mm <- ncol(Z)

  FF <- as.matrix(FF)
  if (nrow(FF) != nn) stop("nrow(FF) must equal the number of cross-section units N = nrow(X)/tt.")
  rr <- ncol(FF)

  if (!is.numeric(nu) || length(nu) != 1L || nu <= 0 || nu >= 1)
    stop("'nu' must be a single number strictly inside (0, 1).")

  ## ---- priors ---------------------------------------------------------------
  fix_prior <- function(pr, nr, nc, name) {
    if (is.null(pr)) return(matrix(1 / nc, nr, nc))
    pr <- as.matrix(pr)
    if (!all(dim(pr) == c(nr, nc))) stop("'", name, "' must be ", nr, " x ", nc, ".")
    if (any(pr <= 0)) stop("All entries of '", name, "' must be strictly positive.")
    rs <- rowSums(pr)
    if (any(abs(rs - 1) > 1e-8)) { warning("Rows of '", name, "' renormalized to sum to 1."); pr <- pr / rs }
    pr
  }
  p0 <- fix_prior(p0, kk, mm, "p0")
  g0 <- fix_prior(g0, nn, rr, "g0")

  ## ---- error support v and its prior w0 -------------------------------------
  ## Default scale: 3 * within-residual sd (Swamy-Arora spirit), not sd(y),
  ## since sd(y) confounds signal variation with the disturbance scale.
  if (is.null(v) || (length(v) == 1L && v == round(v) && v >= 2)) {
    ybar <- ave(y, id); Xbar <- apply(X, 2L, function(col) ave(col, id))
    wres <- tryCatch(residuals(lm.fit(X - Xbar, y - ybar)), error = function(e) y - ybar)
    s_e  <- stats::sd(wres); if (!is.finite(s_e) || s_e == 0) s_e <- stats::sd(y)
    jj <- if (is.null(v)) 3L else as.integer(v)
    v  <- seq(-3 * s_e, 3 * s_e, length.out = jj)
  } else {
    v <- as.numeric(v)
    if (length(v) < 2L) stop("'v' must be an error-support vector (length >= 2) or an integer count of support points.")
    if (abs(mean(range(v))) > 1e-8 * diff(range(v)))
      warning("Error support 'v' is not centered at zero; the implied error prior mean is nonzero.")
  }
  jj <- length(v)
  if (is.null(w0)) {
    w0 <- matrix(1 / jj, nt, jj)
  } else {
    w0 <- as.matrix(w0)
    if (!all(dim(w0) == c(nt, jj))) stop("'w0' must be ", nt, " x ", length(v), " (rows = observations, cols = support points).")
    if (any(w0 <= 0)) stop("All entries of 'w0' must be strictly positive.")
    rs <- rowSums(w0)
    if (any(abs(rs - 1) > 1e-8)) { warning("Rows of 'w0' renormalized to sum to 1."); w0 <- w0 / rs }
  }
  logp0 <- log(p0); logg0 <- log(g0); logw0 <- log(w0)

  ## ---- layout: how lambda (length nt) maps to (n, t) ------------------------
  unit_sum <- function(x) {                    # returns length-nn vector sum_t x_nt
    if (layout == "period") rowSums(matrix(x, nrow = nn, ncol = tt))
    else                    colSums(matrix(x, nrow = tt, ncol = nn))
  }
  unit_rep <- function(m) {                    # replicates length-nn vector to length nt
    if (layout == "period") rep(m, times = tt) else rep(m, each = tt)
  }

  ## ---- stabilized softmax pieces --------------------------------------------
  lse_rows <- function(E, logpr) {             # rows: log sum_j pr_j exp(E_j), and softmax
    L  <- logpr + E
    mx <- apply(L, 1L, max)
    W  <- exp(L - mx); rs <- rowSums(W)
    list(logZ = mx + log(rs), prob = W / rs)
  }

  pieces <- function(lambda) {
    a  <- as.numeric(crossprod(X, lambda))                     # K:  (X'lambda)_k
    sp <- lse_rows(Z * (a / nu), logp0)                        # K x M
    b  <- rowSums(Z * sp$prob)                                 # beta_k

    s  <- unit_sum(lambda)                                     # N:  sum_t lambda_nt
    sg <- lse_rows(FF * (s / nu), logg0)                       # N x R
    mu <- rowSums(FF * sg$prob)                                # mu_n

    sw <- lse_rows(outer(lambda, v) / (1 - nu), logw0)         # NT x J
    e  <- as.numeric(sw$prob %*% v)                            # e_nt

    list(b = b, mu = mu, e = e, p = sp$prob, g = sg$prob, w = sw$prob,
         logOmega = sp$logZ, logPhi = sg$logZ, logPsi = sw$logZ)
  }

  obj <- function(lambda) {
    pc <- pieces(lambda)
    sum(y * lambda) - nu * sum(pc$logOmega) - nu * sum(pc$logPhi) -
      (1 - nu) * sum(pc$logPsi)
  }
  grad <- function(lambda) {                                   # = y - Xb - mu - e
    pc <- pieces(lambda)
    y - as.numeric(X %*% pc$b) - unit_rep(pc$mu) - pc$e
  }

  ## ---- solve -----------------------------------------------------------------
  con <- list(maxit = 1000L, reltol = 1e-12)
  con[names(control)] <- control
  con$fnscale <- -1
  est <- stats::optim(rep(0, nt), fn = obj, gr = grad, method = "BFGS", control = con)
  if (est$convergence != 0L)
    warning("optim did not converge (code ", est$convergence, ").")

  pc <- pieces(est$par)
  fitted <- as.numeric(X %*% pc$b) + unit_rep(pc$mu)
  resid  <- y - fitted

  ent_pos <- function(z) { z <- z[z > 0]; -sum(z * log(z)) }
  Hp <- -sum(pc$p * log(pc$p)); Hg <- -sum(pc$g * log(pc$g)); Hw <- -sum(pc$w * log(pc$w))
  H_signal <- apply(pc$p, 1L, ent_pos)                        # K-vector per-coef entropy
  S_p <- Hp / (kk * log(mm))

  b <- pc$b
  rn <- colnames(X); if (is.null(rn)) rn <- paste0("x", seq_len(kk))
  names(b) <- rn

  ## ---- standard errors (Lee & Cheon 2014, Sec. 3.3) -------------------------
  ## beta:  Var(beta) per vcov_type (within / cluster / disturbance eq. 3.17).
  ## p, mu: delta method from Var(beta). All use only X, residuals, p, beta, so
  ##        they are independent of the fnscale = -1 sign convention.
  cnt    <- as.vector(table(id))                       # periods per unit (T_n)
  xbar_n <- rowsum(X, id) / cnt                        # N x K unit means (units 1..N)
  Xtil   <- X - xbar_n[id, , drop = FALSE]             # within-demeaned design
  eps_res <- resid                                     # idiosyncratic residual
  u_hat   <- y - as.numeric(X %*% b)                   # composite disturbance
  df_e    <- max(nt - nn - kk, 1L)
  sigma2_eps <- sum(eps_res^2) / df_e

  cluster_meat <- function(Xd, r) {                    # sum_n (Xd_n' r_n)(Xd_n' r_n)'
    m <- matrix(0, kk, kk)
    for (g in seq_len(nn)) {
      rows <- which(id == g)
      s <- crossprod(Xd[rows, , drop = FALSE], r[rows])
      m <- m + tcrossprod(s)
    }
    m
  }
  XtXtil <- crossprod(Xtil)
  XtXi_w <- tryCatch(solve(XtXtil), error = function(e) {
    warning("within design (Xtilde'Xtilde) is singular (an intercept or a ",
            "time-invariant regressor?); beta SEs set to NA. Consider ",
            "vcov_type = \"disturbance\".")
    matrix(NA_real_, kk, kk)
  })
  Vbeta <- switch(vcov_type,
    within      = sigma2_eps * XtXi_w,
    cluster     = XtXi_w %*% cluster_meat(Xtil, eps_res) %*% XtXi_w,
    disturbance = {
      XtXi <- tryCatch(solve(crossprod(X)), error = function(e)
        matrix(NA_real_, kk, kk))
      XtXi %*% cluster_meat(X, u_hat) %*% XtXi
    })
  Vbeta <- (Vbeta + t(Vbeta)) / 2                      # symmetrize
  dimnames(Vbeta) <- list(rn, rn)
  se_beta <- sqrt(pmax(diag(Vbeta), 0)); names(se_beta) <- rn

  ## p: SE(p_km) = |p_km (z_km - beta_k)| / Var_p(z_k) * SE(beta_k)
  se_p <- matrix(0, kk, mm, dimnames = list(rn, NULL))
  for (k in seq_len(kk)) {
    s2 <- max(sum(pc$p[k, ] * (Z[k, ] - b[k])^2), 1e-12)
    se_p[k, ] <- abs(pc$p[k, ] * (Z[k, ] - b[k])) / s2 * se_beta[k]
  }

  ## mu: SE(mu_n) = sqrt(sigma2_eps / T_n + xbar_n' Var(beta) xbar_n)
  se_mu <- sqrt(pmax(sigma2_eps / cnt + rowSums((xbar_n %*% Vbeta) * xbar_n), 0))

  structure(
    list(
      # draft names
      b_hat = b, mu_hat = pc$mu, e_hat = pc$e,
      p_hat = pc$p, g_hat = pc$g, w_hat = pc$w, lambda_hat = est$par,
      fitted.values = fitted, residuals = resid,
      foc_residual = max(abs(resid - pc$e)),
      H_p = Hp, H_g = Hg, H_w = Hw, S_p = S_p,
      value = est$value, nu = nu, v = v, layout = layout,
      convergence = est$convergence,
      # canonical aliases
      coefficients = b, H_signal = H_signal, S = S_p,
      objective = est$value, converged = (est$convergence == 0L), method = "dual",
      # standard errors (Sec. 3.3; delta method for p, mu)
      vcov = Vbeta, se_beta = se_beta, se_p = se_p, se_mu = se_mu,
      sigma2_eps = sigma2_eps, vcov_type = vcov_type,
      # stored inputs / dims
      Z = Z, FF = FF, p0 = p0, g0 = g0, w0 = w0, X = X, y = y,
      tt = tt, N = nn, K = kk, M = mm, R = rr, J = jj,
      call = match.call()
    ),
    class = c("panel_gce", "infometrics")
  )
}

## ---- S3 methods ------------------------------------------------------------

#' @export
coef.panel_gce <- function(object, ...) object$coefficients

#' @export
fitted.panel_gce <- function(object, ...) object$fitted.values

#' @export
residuals.panel_gce <- function(object, ...) object$residuals

#' @export
vcov.panel_gce <- function(object, ...) object$vcov

#' @export
print.panel_gce <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Panel GCE (dual, Lee & Cheon 2014), nu =", format(x$nu, digits = digits), "\n")
  cat(strrep("-", 50L), "\n")
  cat(sprintf("  N = %d units   T = %d periods   K = %d   nu = %.3g\n",
              x$N, x$tt, x$K, x$nu))
  cat("beta:", paste(format(x$b_hat, digits = digits), collapse = "  "), "\n")
  cat("se  :", paste(format(x$se_beta, digits = digits), collapse = "  "),
      sprintf("  (%s)\n", x$vcov_type))
  cat(sprintf("mu  : %d effects, range [%.4g, %.4g], sd %.4g\n",
              x$N, min(x$mu_hat), max(x$mu_hat), stats::sd(x$mu_hat)))
  cat(sprintf("S(p) = %s   FOC residual = %s   convergence = %s\n",
              format(x$S, digits = digits), format(x$foc_residual, digits = digits),
              if (x$converged) "yes" else "no"))
  invisible(x)
}

#' @export
summary.panel_gce <- function(object, digits = max(3L, getOption("digits") - 3L),
                              ...) {
  cat("\nPanel GCE (dual, Lee & Cheon 2014, eq. 3.12)\n")
  cat("Call: "); print(object$call); cat("\n")
  cat(sprintf("  N=%d units  T=%d periods  K=%d  M=%d  R=%d  layout=%s  nu=%.3g\n",
              object$N, object$tt, object$K, object$M, object$R, object$layout, object$nu))
  cat("\nResiduals (y - X beta - mu):\n")
  print(summary(object$residuals, digits = digits))
  cat(sprintf("\nCoefficients (beta = Z p; SE type: %s):\n", object$vcov_type))
  z  <- object$b_hat / object$se_beta
  ct <- cbind(Estimate = object$b_hat, `Std. Error` = object$se_beta,
              `z value` = z, `Pr(>|z|)` = 2 * pnorm(-abs(z)))
  printCoefmat(ct, digits = digits, has.Pvalue = TRUE)
  cat(sprintf("\nIndividual effects mu: range [%.4g, %.4g], sd %.4g;  SE(mu) in [%.4g, %.4g]\n",
              min(object$mu_hat), max(object$mu_hat), stats::sd(object$mu_hat),
              min(object$se_mu), max(object$se_mu)))
  cat(sprintf("H(p)=%.4f  H(g)=%.4f  H(w)=%.4f   normalized S(p)=%.4f\n",
              object$H_p, object$H_g, object$H_w, object$S))
  cat(sprintf("max|FOC residual| = %s   convergence: %s\n",
              format(object$foc_residual, digits = digits),
              if (object$converged) "yes (0)" else "no"))
  cat("\nNote: beta SEs use the '", object$vcov_type,
      "' estimator (Lee & Cheon 2014, Sec. 3.3); p and mu SEs\n",
      "are delta-method transforms (mu is a finite-T prediction SE). See ?panel_gce.\n",
      sep = "")
  invisible(object)
}
