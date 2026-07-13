# matrix_ce.R
# Cross-Entropy / Maximum-Entropy matrix balancing (Golan 2008, Section 7.2).
# Recovers an n x m column-stochastic matrix p (each column a distribution over
# the n rows) from column weights x (length m) and row-weighted aggregates
# y (length n), via the concentrated dual in the n row multipliers.
#
# Deliberate standalone estimator (matrix balancing / contingency-table
# reconstruction) -- not a duplicate of any existing function. By user
# instruction it keeps its own argument names (y, x, p0) rather than me()'s.
#
# References:
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145.
#     Section 7.2.

# ---- estimator --------------------------------------------------------------

#' Cross-Entropy Matrix Balancing
#'
#' Recovers a non-negative matrix \eqn{p} (each column a probability
#' distribution) from column weights \code{x} and row-weighted aggregates
#' \code{y} by minimizing the cross-entropy \eqn{D(p\|p0)} relative to a prior
#' \code{p0}, solving the concentrated dual of Golan (2008), Section 7.2. With
#' a uniform \code{p0} this reduces to Maximum Entropy. The model is
#' \deqn{y_i = \sum_j p_{ij} x_j, \qquad \sum_i p_{ij} = 1,}
#' with solution \eqn{\hat p_{ij} = p0_{ij}\exp(\hat\lambda_i x_j)/\Omega_j} and
#' \eqn{\Omega_j(\lambda) = \sum_i p0_{ij}\exp(\lambda_i x_j)}.
#'
#' @details
#' The dual is minimised here (convex) to match the package convention used by
#' \code{\link{inverse_ce}}:
#' \deqn{\ell(\lambda) = -\sum_i \lambda_i y_i + \sum_j \log\Omega_j(\lambda),}
#' minimised over the \eqn{n} row multipliers (equivalent to maximising the
#' column entropies subject to the data). The problem is typically
#' under-determined (\eqn{n\times m} unknowns, \eqn{n + m} constraints), so the
#' returned matrix is the minimum cross-entropy matrix relative to \code{p0}
#' consistent with the data (the maximum-entropy matrix when \code{p0} is
#' uniform) -- not necessarily the matrix that generated \code{y}. To bias the
#' solution toward a known target, pass that target as \code{p0}.
#'
#' The normalized entropy is \eqn{S = H(\hat p) / H(p0)} (Golan §7.5), the ratio
#' of the total Shannon entropy to the prior's entropy. For a uniform \code{p0}
#' (Maximum Entropy) \eqn{H(p0) = m\log n} (the absolute maximum, \eqn{m} columns
#' uniform over \eqn{n} rows), so \eqn{S \in [0, 1]} with \eqn{S = 1} when every
#' column is uniform; for a non-uniform prior \eqn{S} is relative to that prior.
#' \code{summary()} reports the §7.5 information measures built on \eqn{S}: the
#' information index \eqn{I = 1 - S}, pseudo-\eqn{R^2 = 1 - S}, the per-column
#' normalized entropy \eqn{S(p_j) = H(p_j)/H(p0_j)}, and the entropy-ratio test
#' \eqn{W = 2[H(p0) - H(\hat p)] \sim \chi^2_{n-1}} for \eqn{H_0\!: P = P0}.
#'
#' @param y Numeric vector of length n: the row-weighted aggregates
#'   \eqn{y_i = \sum_j p_{ij} x_j}.
#' @param x Numeric vector of length m: the column weights (in \code{print},
#'   the fitted \code{matrix_ce} object).
#' @param p0 Optional n-by-m matrix of strictly-positive prior probabilities
#'   whose columns sum to 1. Defaults to uniform (\code{1/n}).
#' @param control A named list merged over the defaults and passed to
#'   \code{\link[stats]{optim}} (BFGS): \code{maxit} (default 500) and
#'   \code{reltol} (default 1e-12).
#'
#' @return An object of class \code{c("matrix_ce", "infometrics")}: a list with
#'   \code{p}/\code{p_hat} (the n-by-m estimate), \code{lambda}/\code{lambda_hat}
#'   (length-n multipliers), \code{fitted}/\code{fitted.values}
#'   (\eqn{\sum_j p_{ij} x_j}), \code{residuals} (\eqn{y -} fitted),
#'   \code{entropy}/\code{H_signal} (Shannon \eqn{H(\hat p)}), \code{entropy_p0}
#'   (prior entropy \eqn{H(p0)}), \code{S} (prior-relative normalized entropy
#'   \eqn{H(\hat p)/H(p0)}; in [0, 1] for uniform \code{p0}), \code{cross_entropy}
#'   (\eqn{D(\hat p\|p0)}),
#'   \code{objective}/\code{value} (dual objective at the optimum), \code{prior},
#'   \code{converged}, \code{convergence}, \code{method}, \code{y}, \code{x},
#'   \code{n}, \code{m}, and \code{call}.
#'
#' @references Golan, A. (2008). \emph{Information and Entropy Econometrics -
#'   A Review and Synthesis}. Foundations and Trends in Econometrics,
#'   2(1-2), 1-145. Section 7.2.
#'
#' @seealso \code{\link{inverse_ce}} for the formula-interface
#'   pure-moment ME/CE estimator.
#'
#' @examples
#' P <- sweep(matrix(1:6, ncol = 2), 2, colSums(matrix(1:6, ncol = 2)), "/")
#' x <- c(3, 2)
#' y <- as.vector(P %*% x)
#' fit <- matrix_ce(y, x)
#' fit$p
#' summary(fit)
#'
#' @importFrom stats optim pchisq
#' @export
matrix_ce <- function(y, x, p0 = NULL, control = list()) {

  cl <- match.call()
  y <- as.vector(y); x <- as.vector(x)
  n <- length(y); m <- length(x)

  if (!is.numeric(y) || !is.numeric(x)) stop("'y' and 'x' must be numeric.")
  if (anyNA(y) || anyNA(x)) stop("'y' and 'x' must not contain NA.")

  if (is.null(p0)) {
    p0 <- matrix(1 / n, nrow = n, ncol = m)
  } else {
    p0 <- as.matrix(p0)
    if (!all(dim(p0) == c(n, m)))
      stop("dim(p0) must be length(y) x length(x) = ", n, " x ", m, ".")
    if (any(p0 <= 0)) stop("All prior masses 'p0' must be strictly positive.")
    cs <- colSums(p0)
    if (any(abs(cs - 1) > 1e-8)) {
      warning("Columns of 'p0' did not sum to 1; renormalizing.")
      p0 <- sweep(p0, 2, cs, "/")
    }
  }
  logp0 <- log(p0)

  ## ---- stable column partition functions (log-sum-exp) ------------------
  ## A_ij = log p0_ij + lambda_i x_j ; normalize within each column j.
  .pieces <- function(lambda) {
    A  <- logp0 + outer(lambda, x)           # n x m
    mj <- apply(A, 2L, max)                   # length m
    EA <- exp(sweep(A, 2L, mj, "-"))
    cj <- colSums(EA)
    list(logOmega = mj + log(cj), p = sweep(EA, 2L, cj, "/"))
  }
  ## minimised convex dual: ell(lambda) = -sum_i lambda_i y_i + sum_j logOmega_j
  obj  <- function(lambda) -sum(lambda * y) + sum(.pieces(lambda)$logOmega)
  grad <- function(lambda) as.vector(.pieces(lambda)$p %*% x - y)

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
  if (est$convergence != 0L)
    warning("optim did not converge (code ", est$convergence, ").")

  lambda <- est$par
  p <- .pieces(lambda)$p
  dimnames(p) <- list(names(y), names(x))
  names(lambda) <- names(y)
  fitted <- as.vector(p %*% x)

  ## ---- entropy diagnostics (Golan Sec. 7.5) -----------------------------
  ent  <- function(M) -sum(M[M > 0] * log(M[M > 0]))
  Hp   <- ent(p)                                    # Shannon H(p)
  Hp0  <- ent(p0)                                   # prior entropy H(p0)
  cent <- sum(p * log(p / p0))                      # D(p || p0)

  structure(
    list(
      p             = p,
      p_hat         = p,                       # CLAUDE.md canonical alias
      lambda        = lambda,
      lambda_hat    = lambda,                  # CLAUDE.md canonical alias
      fitted        = fitted,
      fitted.values = fitted,                  # lm-consistent alias
      residuals     = as.vector(y - fitted),
      entropy       = Hp,
      H_signal      = Hp,                      # CLAUDE.md canonical field (scalar)
      entropy_p0    = Hp0,                      # prior entropy (for Sec. 7.5 summary)
      S             = if (Hp0 > 0) Hp / Hp0 else NA_real_,   # prior-relative (Sec. 7.5)
      cross_entropy = cent,
      objective     = est$value,               # CLAUDE.md canonical field
      value         = est$value,
      prior         = p0,
      converged     = (est$convergence == 0L), # CLAUDE.md canonical field
      convergence   = est$convergence,
      method        = "dual",                  # CLAUDE.md canonical field
      y             = y,
      x             = x,
      n             = n,
      m             = m,
      call          = cl
    ),
    class = c("matrix_ce", "infometrics")
  )
}

## ----------------------------------------------------------------------------
## S3 methods
## ----------------------------------------------------------------------------

#' @describeIn matrix_ce Estimated probability matrix \eqn{\hat p} (n x m).
#' @param object A \code{matrix_ce} object.
#' @param ... Additional arguments passed to or from methods (currently ignored).
#' @export
coef.matrix_ce <- function(object, ...) object$p

#' @describeIn matrix_ce Fitted row aggregates \eqn{\sum_j p_{ij} x_j}.
#' @export
fitted.matrix_ce <- function(object, ...) object$fitted

#' @describeIn matrix_ce Moment residuals \eqn{y -} fitted.
#' @export
residuals.matrix_ce <- function(object, ...) object$residuals

#' @describeIn matrix_ce Compact printed display.
#' @param digits Number of significant digits to print.
#' @export
print.matrix_ce <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Cross-entropy matrix balancing (", x$n, " x ", x$m, ")\n", sep = "")
  cat("\nEstimated p (columns sum to 1):\n")
  print.default(round(x$p, digits))
  cat("\nH(p) =", format(x$entropy, digits = digits),
      "  normalized S =", format(x$S, digits = digits),
      "  D(p||p0) =", format(x$cross_entropy, digits = digits), "\n")
  cat("max|residual| =", format(max(abs(x$residuals)), digits = digits),
      "  convergence:", x$convergence,
      if (x$convergence == 0L) "(converged)" else "(see ?optim)", "\n")
  invisible(x)
}

#' @describeIn matrix_ce Summary with the Section-7.5 information measures, the
#'   entropy-ratio test for the signal matrix, and a moment-fit table.
#' @export
summary.matrix_ce <- function(object, ...) {
  p <- object$p; p0 <- object$prior
  n <- object$n; m <- object$m
  ent <- function(M) -sum(M[M > 0] * log(M[M > 0]))

  Hp  <- object$entropy; Hp0 <- object$entropy_p0
  S   <- object$S                                # H(p)/H(p0), whole system
  ## per-column normalized entropy  S(p_j) = H(p_j)/H(p0_j)
  Hpj  <- apply(p, 2L, ent); Hp0j <- apply(p0, 2L, ent)
  Sj   <- ifelse(Hp0j > 0, Hpj / Hp0j, NA_real_); names(Sj) <- colnames(p)
  ## entropy-ratio test for H0: P = P0 (lambda = 0):  W = 2[H(p0) - H(p)]
  W   <- 2 * (Hp0 - Hp); dfW <- n - 1L

  mom <- data.frame(Observed = object$y,
                    Fitted   = object$fitted,
                    Residual = object$residuals)
  structure(
    list(call          = object$call,
         dims          = c(n = n, m = m),
         entropy       = Hp,
         entropy_p0    = Hp0,
         cross_entropy = object$cross_entropy,
         S             = S,
         I             = 1 - S,
         pseudo_R2     = 1 - S,
         S_col         = Sj,
         W             = W,
         df            = dfW,
         p_value       = stats::pchisq(W, dfW, lower.tail = FALSE),
         moments       = mom,
         value         = object$value,
         convergence   = object$convergence,
         p             = p),
    class = "summary.matrix_ce"
  )
}

#' @describeIn matrix_ce Print method for the summary object.
#' @export
print.summary.matrix_ce <- function(x, digits = max(3L, getOption("digits") - 3L),
                                     ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Cross-entropy matrix balancing (", x$dims["n"], " x ", x$dims["m"], ")\n",
      sep = "")
  cat("\nEstimated p (columns sum to 1):\n")
  print.default(round(x$p, digits))
  cat("\nMoment fit (y = sum_j p_ij x_j):\n")
  if (nrow(x$moments) > 10L) {
    print(summary(x$moments$Residual, digits = digits))
  } else {
    print(round(x$moments, digits))
  }

  cat("\nInformation measures (Golan 7.5):\n")
  cat("  Normalized entropy S(P) =", format(x$S, digits = digits),
      "   Information index I(P) =", format(x$I, digits = digits),
      "   pseudo-R2 =", format(x$pseudo_R2, digits = digits), "\n")
  cat("  Per-column S(p_j):", paste(format(x$S_col, digits = digits), collapse = "  "), "\n")
  cat(sprintf("  Entropy-ratio test  H0: P = P0 :  W = %s  ~ chi2(%d)  p = %s\n",
              format(x$W, digits = digits), x$df, format(x$p_value, digits = digits)))

  cat("\nH(p) =", format(x$entropy, digits = digits),
      "  H(p0) =", format(x$entropy_p0, digits = digits),
      "  D(p||p0) =", format(x$cross_entropy, digits = digits), "\n")
  cat("Convergence code:", x$convergence,
      if (x$convergence == 0L) "(converged)" else "(see ?optim)", "\n")
  invisible(x)
}
