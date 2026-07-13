# me.R
# Classical Maximum Entropy (ME) and Cross-Entropy (CE) estimators.
# Both the primal (constrained) and dual (concentrated, unconstrained)
# formulations are implemented. The dual is the default solver.
#
# References:
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145.
#     Sections 4.1, 4.2.
#   Jaynes, E.T. (1957). Information theory and statistical mechanics.
#     Physical Review, 106, 620-630.

# ---- ME/CE solver -----------------------------------------------------------

#' Maximum Entropy and Cross-Entropy Estimation
#'
#' Estimates a probability distribution p over K outcomes subject to T moment
#' constraints, using the Maximum Entropy (ME) or Cross-Entropy (CE) principle.
#' The default solver uses the dual concentrated (unconstrained) formulation,
#' which is computationally efficient. The primal (constrained) form is
#' available via \code{method = "primal"} when \pkg{nloptr} is installed.
#'
#' @param y Numeric vector of length T. The observed moment values
#'   (right-hand sides of the constraints).
#' @param X Numeric matrix of dimension T x K (or a vector of length K for
#'   T = 1). The moment matrix: column k contains the values x_{1k}, ..., x_{Tk}
#'   such that the constraints are sum_k p_k * x_{tk} = y_t, t = 1,...,T.
#' @param q Numeric vector of length K. Prior (reference) probabilities for
#'   CE estimation. Defaults to the uniform distribution (ME formulation).
#' @param method Character. Either \code{"dual"} (default, concentrated
#'   unconstrained model solved via BFGS) or \code{"primal"} (constrained
#'   optimization via \pkg{nloptr}, must be installed).
#' @param control A named list of control parameters passed to the optimizer:
#'   \describe{
#'     \item{\code{maxit}}{Maximum iterations (default 500).}
#'     \item{\code{tol}}{Convergence tolerance (default 1e-10).}
#'     \item{\code{trace}}{Logical; print optimization progress (default FALSE).}
#'   }
#'
#' @return An object of class \code{"infometrics_me"} (which inherits from
#'   \code{"infometrics"}), a list containing:
#'   \describe{
#'     \item{\code{p_hat}}{Numeric vector of length K: estimated probabilities.}
#'     \item{\code{lambda_hat}}{Numeric vector of length T: optimal Lagrange
#'       multipliers (dual variables). The marginal information of each moment.}
#'     \item{\code{H_star}}{Numeric. Entropy at the solution H(p_hat).}
#'     \item{\code{S}}{Numeric. Normalized entropy S(p_hat) in [0, 1].}
#'     \item{\code{objective}}{Numeric. Value of the objective at convergence.}
#'     \item{\code{converged}}{Logical. Whether the optimizer converged.}
#'     \item{\code{method}}{Character. Solver used ("dual" or "primal").}
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{y}}{The moment vector supplied.}
#'     \item{\code{X}}{The moment matrix supplied.}
#'     \item{\code{q}}{The prior distribution used.}
#'   }
#'
#' @details
#' \strong{ME formulation (uniform priors):}
#' Maximizes H(p) = -sum(p_k * log(p_k)) subject to X'p = y and sum(p_k) = 1.
#'
#' \strong{CE formulation (non-uniform priors):}
#' Minimizes D(p||q) = sum(p_k * log(p_k / q_k)) subject to the same constraints.
#' The ME formulation is a special case with q_k = 1/K for all k.
#'
#' \strong{Dual (concentrated) model:}
#' Both problems reduce to an unconstrained minimization over the T-dimensional
#' Lagrange multiplier vector lambda (Golan, 2008, Eq. 4.4--4.5):
#'
#'   min_lambda  { -lambda'y + log(sum_k q_k * exp(sum_t lambda_t * x_{tk})) }
#'
#' The solution lambda_hat yields the estimated probabilities via:
#'   p_hat_k = q_k * exp(sum_t lambda_hat_t * x_{tk}) / Omega(lambda_hat)
#'
#' where Omega(lambda) = sum_k q_k * exp(sum_t lambda_t * x_{tk}) is the
#' partition function. The dual reduces problem dimensionality from K >> 1
#' to T, the number of moment constraints.
#'
#' @references
#' Golan, A. (2008). Information and Entropy Econometrics. \emph{Foundations
#' and Trends in Econometrics}, \strong{2}(1-2), 1-145.
#'
#' Jaynes, E.T. (1957). Information theory and statistical mechanics.
#' \emph{Physical Review}, \strong{106}, 620-630.
#'
#' @examples
#' # Classic "loaded die" example (Golan 2008, Section 4.6):
#' # We observe only the mean of a 6-sided die and want to recover p.
#' support <- 1:6
#' y_mean  <- 4.5              # observed mean (slightly above fair die mean of 3.5)
#' X       <- matrix(support, nrow = 1)   # 1 moment x 6 outcomes
#'
#' fit <- me(y = y_mean, X = X)
#' print(fit)
#' summary(fit)
#'
#' # With a non-uniform prior (CE formulation)
#' q_prior <- c(0.10, 0.10, 0.15, 0.20, 0.25, 0.20)
#' fit_ce <- me(y = y_mean, X = X, q = q_prior)
#' print(fit_ce)
#'
#' @export
me <- function(y, X, q = NULL, method = c("dual", "primal"), control = list()) {
  mc <- match.call()
  method <- match.arg(method)

  # -- coerce and validate inputs -------------------------------------------
  y <- as.numeric(y)
  T_mom <- length(y)  # number of moment constraints

  X <- as.matrix(X)
  if (nrow(X) != T_mom)
    stop("nrow(X) must equal length(y). Got nrow(X) = ", nrow(X),
         ", length(y) = ", T_mom, ".")
  K <- ncol(X)  # number of support points / outcomes

  if (K <= T_mom)
    stop("The problem must be under-determined: need K > T + 1. ",
         "Got K = ", K, ", T = ", T_mom, ".")

  # default prior: uniform
  if (is.null(q)) {
    q <- rep(1 / K, K)
  } else {
    q <- as.numeric(q)
    .check_prob(q, name = "q")
    if (length(q) != K)
      stop("q must have length K = ", K, ". Got length(q) = ", length(q), ".")
  }

  # -- control defaults -----------------------------------------------------
  con <- list(maxit = 500L, tol = 1e-10, trace = FALSE)
  con[names(control)] <- control

  # -- dispatch to solver ---------------------------------------------------
  if (method == "dual") {
    out <- .me_dual(y, X, q, con)
  } else {
    out <- .me_primal(y, X, q, con)
  }

  out$call <- mc
  out$y    <- y
  out$X    <- X
  out$q    <- q
  class(out) <- c("infometrics_me", "infometrics")
  out
}

# ---- dual (concentrated) solver --------------------------------------------

#' @keywords internal
.me_dual <- function(y, X, q, con) {
  T_mom <- length(y)
  K     <- ncol(X)

  # Partition function: Omega(lambda) = sum_k q_k * exp(X' %*% lambda)
  # Exponents: e_k(lambda) = sum_t lambda_t * x_{tk}  =>  X %*% lambda (K-vec)
  # (X is T x K, lambda is T x 1, so X %*% lambda... but we want K-vector)
  # Using t(X) %*% lambda  [K x T] %*% [T x 1] = K x 1
  .log_omega <- function(lambda) {
    exponents <- as.numeric(t(X) %*% lambda)  # K-vector
    # numerically stable log-sum-exp:
    max_e <- max(exponents)
    log(sum(q * exp(exponents - max_e))) + max_e
  }

  # Concentrated (dual) objective: minimize  -lambda'y + log(Omega(lambda))
  # Equivalent to maximizing the entropy (Golan 2008, Eq. 4.5)
  .dual_obj <- function(lambda) {
    -as.numeric(crossprod(lambda, y)) + .log_omega(lambda)
  }

  # Analytic gradient: d/d(lambda_t) = -y_t + sum_k p_k(lambda) * x_{tk}
  # i.e. gradient = -(y - X %*% p(lambda))  [residual of moment constraints]
  .dual_grad <- function(lambda) {
    exponents <- as.numeric(t(X) %*% lambda)
    max_e     <- max(exponents)
    exp_adj   <- q * exp(exponents - max_e)
    p_lam     <- exp_adj / sum(exp_adj)        # K-vector of probabilities
    -y + as.numeric(X %*% p_lam)              # T-vector gradient
  }

  # Optimize via BFGS
  lambda_init <- rep(0, T_mom)
  opt <- stats::optim(
    par     = lambda_init,
    fn      = .dual_obj,
    gr      = .dual_grad,
    method  = "BFGS",
    control = list(
      maxit   = con$maxit,
      reltol  = con$tol,
      trace   = as.integer(con$trace)
    )
  )

  # Recover probabilities from optimal lambda
  exponents <- as.numeric(t(X) %*% opt$par)
  max_e     <- max(exponents)
  exp_adj   <- q * exp(exponents - max_e)
  p_hat     <- exp_adj / sum(exp_adj)

  H_star <- shannon_entropy(p_hat)
  S      <- normalized_entropy(p_hat, q)

  list(
    p_hat      = p_hat,
    lambda_hat = opt$par,
    H_star     = H_star,
    S          = S,
    objective  = opt$value,
    converged  = (opt$convergence == 0),
    method     = "dual"
  )
}

# ---- primal (constrained) solver -------------------------------------------

#' @keywords internal
.me_primal <- function(y, X, q, con) {
  if (!requireNamespace("nloptr", quietly = TRUE))
    stop("Package 'nloptr' is required for method = 'primal'. ",
         "Install it with: install.packages('nloptr')")

  K     <- ncol(X)
  T_mom <- length(y)

  # Objective: minimize KL divergence D(p||q) = sum p_k log(p_k / q_k)
  # (Maximizing H(p) for ME is equivalent to minimizing D(p||uniform))
  .obj <- function(p) {
    idx <- p > 0
    sum(p[idx] * log(p[idx] / q[idx]))
  }
  .obj_grad <- function(p) {
    g <- log(p / q) + 1
    g[p <= 0] <- -Inf
    g
  }

  # Equality constraints: X %*% p - y = 0  and  sum(p) - 1 = 0
  .heq <- function(p) c(as.numeric(X %*% p) - y, sum(p) - 1)
  .heq_jac <- function(p) rbind(X, rep(1, K))

  # Starting point: uniform (or prior)
  p_init <- q

  res <- nloptr::nloptr(
    x0          = p_init,
    eval_f      = .obj,
    eval_grad_f = .obj_grad,
    lb          = rep(0, K),
    ub          = rep(1, K),
    eval_g_eq   = .heq,
    eval_jac_g_eq = .heq_jac,
    opts = list(
      algorithm   = "NLOPT_LD_SLSQP",
      maxeval     = con$maxit,
      xtol_rel    = con$tol,
      print_level = if (con$trace) 1L else 0L
    )
  )

  p_hat  <- res$solution
  p_hat  <- pmax(p_hat, 0)
  p_hat  <- p_hat / sum(p_hat)  # renormalize for numerical safety

  # Recover Lagrange multipliers via the primal-dual relationship
  # (Approximate: solve for lambda from the FOCs)
  lambda_hat <- tryCatch({
    as.numeric(
      solve(X %*% t(X), X %*% (log(p_hat / q) + 1))
    )
  }, error = function(e) rep(NA_real_, T_mom))

  H_star <- shannon_entropy(p_hat)
  S      <- normalized_entropy(p_hat, q)

  list(
    p_hat      = p_hat,
    lambda_hat = lambda_hat,
    H_star     = H_star,
    S          = S,
    objective  = res$objective,
    converged  = (res$status > 0),
    method     = "primal"
  )
}

# ---- S3 methods -------------------------------------------------------------

#' Print method for infometrics_me objects
#'
#' @param x An object of class \code{"infometrics_me"}.
#' @param digits Integer. Number of significant digits (default 4).
#' @param ... Not used.
#' @export
print.infometrics_me <- function(x, digits = 4, ...) {
  cat("Maximum Entropy / Cross-Entropy Estimator\n")
  cat(rep("-", 45), "\n", sep = "")
  cat("Method      :", x$method, "\n")
  cat("Outcomes  K :", length(x$p_hat), "\n")
  cat("Moments   T :", length(x$y), "\n")
  cat("Converged   :", x$converged, "\n")
  cat(rep("-", 45), "\n", sep = "")
  cat("Estimated probabilities:\n")
  print(round(x$p_hat, digits))
  cat("\nLagrange multipliers:\n")
  print(round(x$lambda_hat, digits))
  cat(rep("-", 45), "\n", sep = "")
  cat("H*(p_hat)           :", round(x$H_star, digits), "\n")
  cat("S (normalized H)    :", round(x$S,      digits), "\n")
  cat("Pseudo-R2 (1 - S)   :", round(1 - x$S,  digits), "\n")
  invisible(x)
}

#' Summary method for infometrics_me objects
#'
#' @param object An object of class \code{"infometrics_me"}.
#' @param digits Integer. Number of significant digits (default 4).
#' @param ... Not used.
#' @export
summary.infometrics_me <- function(object, digits = 4, ...) {
  cat("=== Maximum Entropy / Cross-Entropy Estimation ===\n\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")

  K <- length(object$p_hat)
  T_mom <- length(object$y)
  cat("Dimensions: K =", K, "outcomes,  T =", T_mom, "moment constraints\n")
  cat("Solver    :", object$method, "\n")
  cat("Converged :", object$converged, "\n\n")

  # Moment fit table
  y_hat <- as.numeric(object$X %*% object$p_hat)
  resid <- object$y - y_hat
  mom_tab <- data.frame(
    Observed  = round(object$y,   digits),
    Fitted    = round(y_hat,      digits),
    Residual  = round(resid,      digits)
  )
  rownames(mom_tab) <- paste0("Moment ", seq_len(T_mom))
  cat("Moment constraints:\n")
  print(mom_tab)
  cat("\n")

  # Entropy-ratio test statistic (chi-squared)
  # W = 2 * K * log(K) * (1 - S)  under uniform priors
  # Under H0: W ~ chi^2(K-1)
  W   <- 2 * K * log(K) * (1 - object$S)
  df  <- K - 1
  pval <- stats::pchisq(W, df = df, lower.tail = FALSE)

  cat("Goodness of fit:\n")
  cat("  H*(p_hat)             =", round(object$H_star, digits), "\n")
  cat("  H_max  = log(K)       =", round(log(K),        digits), "\n")
  cat("  S (normalized entropy)=", round(object$S,      digits), "\n")
  cat("  Pseudo-R2  (1 - S)    =", round(1 - object$S,  digits), "\n\n")
  cat("Entropy-ratio test (H0: all parameters = 0):\n")
  cat("  W  =", round(W, digits), "\n")
  cat("  df =", df, "\n")
  cat("  p-value =", format.pval(pval, digits = digits), "\n\n")

  # Lagrange multipliers
  cat("Lagrange multipliers (marginal information of each moment):\n")
  lam_tab <- data.frame(
    lambda = round(object$lambda_hat, digits)
  )
  rownames(lam_tab) <- paste0("lambda[", seq_len(T_mom), "]")
  print(lam_tab)

  invisible(object)
}

#' Extract estimated probabilities from an infometrics_me object
#'
#' @param object An object of class \code{"infometrics_me"}.
#' @param ... Not used.
#' @return Named numeric vector of estimated probabilities.
#' @export
coef.infometrics_me <- function(object, ...) {
  p <- object$p_hat
  names(p) <- paste0("p[", seq_along(p), "]")
  p
}

#' Fitted (moment) values from an infometrics_me object
#'
#' @param object An object of class \code{"infometrics_me"}.
#' @param ... Not used.
#' @return Numeric vector of fitted moment values X %*% p_hat.
#' @export
fitted.infometrics_me <- function(object, ...) {
  as.numeric(object$X %*% object$p_hat)
}

#' Residuals from an infometrics_me object
#'
#' Returns the moment residuals y - X %*% p_hat. For the dual solver these
#' should be at machine-precision zero at convergence.
#'
#' @param object An object of class \code{"infometrics_me"}.
#' @param ... Not used.
#' @return Numeric vector of moment residuals.
#' @export
residuals.infometrics_me <- function(object, ...) {
  object$y - stats::fitted(object)
}
