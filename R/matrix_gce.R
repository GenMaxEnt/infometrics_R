# matrix_gce.R
# Generalized Cross-Entropy matrix balancing with stochastic moments
# (Golan 2008, Section 7.4). Noisy extension of matrix_ce(): the aggregates
# are y_i = sum_j p_ij x_j + e_i, with signal p (columns sum to 1) and noise
# w (rows sum to 1) recovered jointly via the concentrated dual.
#
# Deliberate standalone estimator (the noisy sibling of matrix_ce, just as
# inverse_noise is the noisy sibling of inverse_ce). By user instruction it
# keeps its own argument names (y, x, p0, w0, v, nu) rather than me()'s.
#
# References:
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145.
#     Section 7.4 (stochastic moments), Section 7.5 (information measures).

# ---- estimator --------------------------------------------------------------

#' Generalized Cross-Entropy Matrix Balancing (Stochastic Moments)
#'
#' Noisy extension of \code{\link{matrix_ce}} (Golan 2008, Section 7.4). The
#' aggregates are treated as stochastic, \eqn{y_i = \sum_j p_{ij} x_j + e_i}
#' with \eqn{e_i = \sum_h v_h w_{ih}}, and the cross-entropy of both the signal
#' \code{p} (columns sum to 1) and the noise \code{w} (rows sum to 1) is
#' minimized relative to priors \code{p0}, \code{w0}.
#'
#' @details
#' The dual is minimised here (convex) to match the package convention used by
#' \code{\link{gme}}, \code{\link{inverse_noise}} and \code{\link{matrix_ce}}:
#' \deqn{
#'   \ell(\lambda) = -\sum_i \lambda_i y_i + \nu \sum_j \log\Omega_j(\lambda)
#'        + (1-\nu)\sum_i \log\Psi_i(\lambda),
#' }
#' minimised over the \eqn{n} row multipliers, with
#' \eqn{\hat p_{ij} \propto p0_{ij}\exp(\lambda_i x_j/\nu)} (normalised within
#' each column) and \eqn{\hat w_{ih} \propto w0_{ih}\exp(\lambda_i v_h/(1-\nu))}
#' (normalised within each row); \eqn{\hat e_i = \sum_h v_h \hat w_{ih}}.
#' Equivalent to maximising the joint signal-plus-noise entropy subject to the
#' stochastic-moment constraints.
#'
#' The weight \eqn{\nu \in (0,1)} (an extension to Eq. 7.22) balances prediction
#' against precision by trading signal entropy against noise entropy. At
#' \eqn{\nu = 0.5} this is the equal-weight GCE of Eq. (7.22) and the estimates
#' coincide with it (only the \eqn{\lambda} scale differs). As the noise support
#' \code{v} shrinks toward 0 the noise is forced out and the solution returns to
#' \code{\link{matrix_ce}}.
#'
#' @param y,x Numeric vectors (lengths n and m), each normalized to [0, 1] so
#'   the implied errors lie in [-1, 1].
#' @param p0 Optional n-by-m signal prior, columns summing to 1, strictly
#'   positive. Default uniform.
#' @param w0 Optional n-by-H noise prior, rows summing to 1, strictly positive.
#'   Default uniform; its column count must match \code{length(v)}.
#' @param v Noise support, symmetric around 0 in [-1, 1]: \code{NULL} (default
#'   \code{c(-1, 0, 1)}), a single whole number read as a support-point count
#'   \code{H >= 2}, or an explicit numeric vector (used as given).
#' @param nu Entropy weight in \eqn{(0, 1)}; default 0.5.
#' @param control A named list merged over the defaults and passed to
#'   \code{\link[stats]{optim}} (BFGS): \code{maxit} (default 500) and
#'   \code{reltol} (default 1e-12).
#'
#' @return An object of class \code{c("matrix_gce", "infometrics")}: a list with
#'   \code{p}/\code{p_hat} (n-by-m signal), \code{w}/\code{w_hat} (n-by-H noise),
#'   \code{e} (length-n noise), \code{lambda}/\code{lambda_hat}, \code{nu},
#'   \code{fitted.values} (\eqn{Px}), \code{residuals} (\eqn{y - Px}, equal to
#'   \code{e} at the optimum), \code{entropy_p}/\code{H_signal} (scalar
#'   \eqn{H(\hat p)}), \code{entropy_w}, \code{entropy_p0}, \code{entropy_w0},
#'   \code{S_p}/\code{S} and \code{S_w} (normalized entropies
#'   \eqn{H(\cdot)/H(\cdot_0)} in [0, 1]), \code{objective}/\code{value} (dual
#'   objective), \code{signal_prior}, \code{noise_prior}, \code{support},
#'   \code{converged}, \code{convergence}, \code{method}, and \code{call}.
#'
#' @references Golan, A. (2008). \emph{Information and Entropy Econometrics -
#'   A Review and Synthesis}. Foundations and Trends in Econometrics,
#'   2(1-2), 1-145. Section 7.4.
#'
#' @seealso \code{\link{matrix_ce}} (the exact-moment version this reduces to as
#'   \code{v} shrinks to 0).
#'
#' @examples
#' set.seed(1)
#' P <- matrix(runif(6), ncol = 2); P <- sweep(P, 2, colSums(P), "/")
#' x <- c(0.6, 0.3); y <- as.vector(P %*% x)
#' fit <- matrix_gce(y, x, nu = 0.5)
#' coef(fit)        # estimated signal matrix p
#' residuals(fit)   # y - Px  (= estimated noise e)
#' summary(fit)
#'
#' @importFrom stats optim pchisq
#' @export
matrix_gce <- function(y, x, p0 = NULL, w0 = NULL, v = NULL, nu = 0.5,
                       control = list()) {

  cl <- match.call()
  y <- as.vector(y); x <- as.vector(x)
  n <- length(y); m <- length(x)

  if (!is.numeric(y) || !is.numeric(x)) stop("'y' and 'x' must be numeric.")
  if (anyNA(y) || anyNA(x)) stop("'y' and 'x' must not contain NA.")
  if (min(y) < 0 || max(y) > 1 || min(x) < 0 || max(x) > 1)
    warning("'y' and 'x' should be normalized to [0, 1] (Golan 7.4); the noise ",
            "support assumes errors in [-1, 1].")
  if (length(nu) != 1L || !is.finite(nu) || nu <= 0 || nu >= 1)
    stop("'nu' must be a single value strictly between 0 and 1.")

  ## ---- signal prior p0 (n x m, columns sum to 1) ------------------------
  if (is.null(p0)) {
    p0 <- matrix(1 / n, n, m)
  } else {
    p0 <- as.matrix(p0)
    if (!all(dim(p0) == c(n, m)))
      stop("dim(p0) must be length(y) x length(x) = ", n, " x ", m, ".")
    if (any(p0 <= 0)) stop("All signal prior masses 'p0' must be positive.")
    cs <- colSums(p0)
    if (any(abs(cs - 1) > 1e-8)) { warning("Columns of 'p0' renormalized."); p0 <- sweep(p0, 2, cs, "/") }
  }

  ## ---- noise support v and prior w0 (n x H, rows sum to 1) --------------
  is_count <- function(z) length(z) == 1L && is.finite(z) && abs(z - round(z)) < 1e-8
  if (is.null(v)) {
    H <- if (is.null(w0)) 3L else ncol(w0)
    v <- seq(-1, 1, length.out = H)
  } else if (is_count(v)) {
    H <- as.integer(round(v)); if (H < 2L) stop("Scalar 'v' is a support count and must be >= 2.")
    v <- seq(-1, 1, length.out = H)
  } else {
    v <- as.vector(v); H <- length(v); if (H < 2L) stop("'v' needs >= 2 points.")
  }
  if (abs(mean(v)) > 1e-8) warning("Noise support 'v' is not symmetric around 0.")
  if (min(v) < -1 - 1e-8 || max(v) > 1 + 1e-8)
    warning("Noise support 'v' extends beyond [-1, 1].")
  if (is.null(w0)) {
    w0 <- matrix(1 / H, n, H)
  } else {
    w0 <- as.matrix(w0)
    if (ncol(w0) != H) stop("ncol(w0) (", ncol(w0), ") must equal length(v) (", H, ").")
    if (nrow(w0) != n) stop("nrow(w0) (", nrow(w0), ") must equal length(y) (", n, ").")
    if (any(w0 <= 0)) stop("All noise prior masses 'w0' must be positive.")
    rw <- rowSums(w0)
    if (any(abs(rw - 1) > 1e-8)) { warning("Rows of 'w0' renormalized."); w0 <- sweep(w0, 1, rw, "/") }
  }

  logp0 <- log(p0); logw0 <- log(w0)

  ## ---- stable pieces (log-sum-exp; signal by column, noise by row) ------
  .signal <- function(lambda) {
    A  <- logp0 + outer(lambda, x) / nu              # n x m
    mj <- apply(A, 2L, max); EA <- exp(sweep(A, 2L, mj, "-")); cj <- colSums(EA)
    list(logOmega = mj + log(cj), p = sweep(EA, 2L, cj, "/"))
  }
  .noise <- function(lambda) {
    B  <- logw0 + outer(lambda, v) / (1 - nu)        # n x H
    mi <- apply(B, 1L, max); EB <- exp(sweep(B, 1L, mi, "-")); ri <- rowSums(EB)
    w  <- sweep(EB, 1L, ri, "/")
    list(logPsi = mi + log(ri), w = w, e = as.vector(w %*% v))
  }
  ## minimised convex dual: -sum lambda_i y_i + nu sum logOmega_j + (1-nu) sum logPsi_i
  obj <- function(lambda) {
    sg <- .signal(lambda); ns <- .noise(lambda)
    -sum(lambda * y) + nu * sum(sg$logOmega) + (1 - nu) * sum(ns$logPsi)
  }
  grad <- function(lambda) {
    sg <- .signal(lambda); ns <- .noise(lambda)
    as.vector(sg$p %*% x + ns$e - y)
  }

  ## ---- minimise the dual -------------------------------------------------
  con <- list(maxit = 500L, reltol = 1e-12)
  con[names(control)] <- control
  est <- stats::optim(
    par     = rep(0, n),
    fn      = obj,
    gr      = grad,
    method  = "BFGS",
    control = list(maxit = con$maxit, reltol = con$reltol)
  )
  if (est$convergence != 0L) warning("optim did not converge (code ", est$convergence, ").")

  lambda <- est$par; names(lambda) <- names(y)
  sg <- .signal(lambda); ns <- .noise(lambda)
  p <- sg$p; dimnames(p) <- list(names(y), names(x))
  w <- ns$w; e <- ns$e
  fitted <- as.vector(p %*% x); names(fitted) <- names(y)

  ent <- function(M) -sum(M[M > 0] * log(M[M > 0]))
  Hp <- ent(p); Hw <- ent(w); Hp0 <- ent(p0); Hw0 <- ent(w0)
  S_p <- if (Hp0 > 0) Hp / Hp0 else NA_real_
  S_w <- if (Hw0 > 0) Hw / Hw0 else NA_real_

  structure(
    list(
      p             = p,
      p_hat         = p,                       # CLAUDE.md canonical alias
      w             = w,
      w_hat         = w,                       # CLAUDE.md canonical alias
      e             = e,
      lambda        = lambda,
      lambda_hat    = lambda,                  # CLAUDE.md canonical alias
      nu            = nu,
      fitted.values = fitted,
      residuals     = as.vector(y - fitted),
      entropy_p     = Hp,
      H_signal      = Hp,                      # CLAUDE.md canonical field (scalar)
      entropy_w     = Hw,
      entropy_p0    = Hp0,
      entropy_w0    = Hw0,
      S_p           = S_p,
      S             = S_p,                     # CLAUDE.md canonical field (signal norm. entropy)
      S_w           = S_w,
      objective     = est$value,               # CLAUDE.md canonical field
      value         = est$value,
      signal_prior  = p0,
      noise_prior   = w0,
      support       = v,
      converged     = (est$convergence == 0L), # CLAUDE.md canonical field
      convergence   = est$convergence,
      method        = "dual",                  # CLAUDE.md canonical field
      call          = cl
    ),
    class = c("matrix_gce", "infometrics")
  )
}

## ----------------------------------------------------------------------------
## S3 methods
## ----------------------------------------------------------------------------

#' @describeIn matrix_gce Estimated signal matrix \code{p}.
#' @param object,x A \code{matrix_gce} object (\code{x} in \code{print}).
#' @export
coef.matrix_gce <- function(object, ...) object$p

#' @describeIn matrix_gce Fitted aggregates \eqn{Px}.
#' @export
fitted.matrix_gce <- function(object, ...) object$fitted.values

#' @describeIn matrix_gce Residuals \eqn{y - Px} (the estimated noise \code{e}).
#' @export
residuals.matrix_gce <- function(object, ...) object$residuals

#' @describeIn matrix_gce Compact display.
#' @param digits Significant digits to print.
#' @export
print.matrix_gce <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("GCE matrix balancing, stochastic moments (", nrow(x$p), " x ",
      ncol(x$p), "),  nu = ", format(x$nu, digits = digits), "\n", sep = "")
  cat("\nCall: ", paste(deparse(x$call), collapse = " "), "\n", sep = "")
  cat("\nEstimated signal p (columns sum to 1):\n"); print.default(round(x$p, digits))
  cat("\nH(p) =", format(x$entropy_p, digits = digits),
      "  H(w) =", format(x$entropy_w, digits = digits),
      "  max|y - Px - e| =", format(max(abs(x$residuals - x$e)), digits = digits), "\n")
  cat("convergence:", x$convergence,
      if (x$convergence == 0L) "(converged)" else "(see ?optim)", "\n")
  invisible(x)
}

#' @describeIn matrix_gce Summary with the Section-7.5 information measures
#'   and entropy-ratio test for the signal matrix.
#' @export
summary.matrix_gce <- function(object, ...) {
  p <- object$p; p0 <- object$signal_prior
  n <- nrow(p); m <- ncol(p)
  ent <- function(M) -sum(M[M > 0] * log(M[M > 0]))

  Hp  <- object$entropy_p; Hp0 <- object$entropy_p0
  S   <- object$S_p                              # H(p)/H(p0), whole system
  ## per-column normalized entropy  S(p_j) = H(p_j)/H(p0_j)
  Hpj  <- apply(p, 2L, ent); Hp0j <- apply(p0, 2L, ent)
  Sj   <- ifelse(Hp0j > 0, Hpj / Hp0j, NA_real_); names(Sj) <- colnames(p)
  ## entropy-ratio test for H0: P = P0 (lambda = 0):  W = 2[H(p0) - H(p)]
  W   <- 2 * (Hp0 - Hp); dfW <- n - 1L

  structure(
    list(call = object$call,
         dims = c(n = n, m = m, H = length(object$support)),
         nu = object$nu, p = p, e = object$e,
         entropy_p = Hp, entropy_w = object$entropy_w,
         S = S, I = 1 - S, pseudo_R2 = 1 - S, S_col = Sj,
         W = W, df = dfW, p_value = stats::pchisq(W, dfW, lower.tail = FALSE),
         value = object$value,
         constraint = max(abs(object$residuals - object$e)),
         convergence = object$convergence),
    class = "summary.matrix_gce")
}

#' @describeIn matrix_gce Print method for the summary object.
#' @export
print.summary.matrix_gce <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat(sprintf("GCE matrix balancing: %d x %d signal, %d noise points, nu = %s\n",
              x$dims["n"], x$dims["m"], x$dims["H"], format(x$nu, digits = digits)))
  cat("\nSignal p (columns sum to 1):\n"); print.default(round(x$p, digits))
  cat("\nNoise e (row errors):\n"); print(summary(x$e, digits = digits))

  cat("\nInformation measures (Golan 7.5):\n")
  cat("  Normalized entropy S(P) =", format(x$S, digits = digits),
      "   Information index I(P) =", format(x$I, digits = digits),
      "   pseudo-R2 =", format(x$pseudo_R2, digits = digits), "\n")
  cat("  Per-column S(p_j):", paste(format(x$S_col, digits = digits), collapse = "  "), "\n")
  cat(sprintf("  Entropy-ratio test  H0: P = P0 :  W = %s  ~ chi2(%d)  p = %s\n",
              format(x$W, digits = digits), x$df, format(x$p_value, digits = digits)))

  cat("\nDual objective:", format(x$value, digits = digits),
      "  H(p) =", format(x$entropy_p, digits = digits),
      "  H(w) =", format(x$entropy_w, digits = digits),
      "  max|y - Px - e|:", format(x$constraint, digits = digits), "\n")
  cat("Convergence code:", x$convergence,
      if (x$convergence == 0L) "(converged)" else "(see ?optim)", "\n")
  invisible(x)
}
