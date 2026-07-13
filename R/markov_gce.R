# markov_gce.R
# Generalized Cross-Entropy (GCE) estimation of a first-order Markov
# transition / balancing matrix from a balanced long panel, with stochastic
# (noisy) moments -- Golan (2008), Section 7.7.1, eqs. (7.25a) and (7.28).
#
# Each conditional moment carries an additive error eps_itj = sum_m w_itjm v_m
# on a support v symmetric on [-1, 1]; this makes the (generally over-
# identified) time-varying-covariate moments feasible. The nu-weighted
# concentrated dual is maximised (fnscale = -1):
#   l(lambda) = sum_{t>=2} sum_j sum_{i,s} y_itj z_its lambda_sj
#               - nu      * sum_k       log Omega_k(lambda)
#               - (1 - nu)* sum_{i,t,j} log Psi_itj(lambda)
# As nu -> 1 / v -> 0 the noise is suppressed and markov_gce -> markov_ce.
#
# Noisy sibling of markov_ce (cf. matrix_gce vs matrix_ce). Reuses the panel
# reshaper .markov_arrays() from markov_ce.R. By user instruction it keeps its
# own argument names (data, id, time, states, covariates, v, nu, p0, w0). The
# dual is kept in the maximise form (fnscale = -1) by user choice -- one of the
# two maximise-form estimators in the package (with markov_ce).
#
# References:
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145.
#     Section 7.7.1, eqs. (7.25a), (7.28).

#' Generalized Cross-Entropy Estimation of a Markov Matrix with Noisy Moments
#'
#' Estimates a row-stochastic \eqn{K \times K} matrix \eqn{P} linking simplex-
#' valued state vectors over time, conditioning on covariates, by Generalized
#' Cross-Entropy from a balanced long panel (Golan 2008, Section 7.7.1, eqs.
#' 7.25a and 7.28). Unlike \code{\link{markov_ce}} (which imposes the moments
#' exactly), each conditional moment carries an additive error
#' \eqn{\varepsilon_{itj} = \sum_m w_{itjm} v_m} on a support \eqn{v} symmetric
#' on \eqn{[-1, 1]}. This makes over-identified moment systems -- in particular
#' time-varying covariates such as household income -- feasible, where the exact
#' problem is unbounded.
#'
#' @details
#' The \eqn{\nu}-weighted concentrated dual maximized over the \eqn{S \times K}
#' multipliers \eqn{\lambda} is
#' \deqn{\ell(\lambda) = \sum_{t\ge 2}\sum_j\sum_{i,s} y_{itj} z_{its}\lambda_{sj}
#'   - \nu \sum_k \log\Omega_k(\lambda)
#'   - (1-\nu)\sum_{i,t,j}\log\Psi_{itj}(\lambda),}
#' with \eqn{\hat p_{kj}\propto p0_{kj}\exp(\nu^{-1}\sum A_{ks}\lambda_{sj})} and
#' \eqn{\hat w_{itjm}\propto w0_{m}\exp((1-\nu)^{-1}\sum_s z_{its}v_m\lambda_{sj})}.
#' The first-order conditions are the noisy moments \eqn{M - A'\hat P -
#' \sum_{i,t} z\,\hat\varepsilon = 0}. The exact estimator \code{\link{markov_ce}}
#' is recovered in two limits: as the support \eqn{v \to 0} (the slack vanishes
#' mechanically), and as \eqn{\nu \to 0} (the noise-divergence term dominates,
#' forcing \eqn{\hat\varepsilon \to 0}). Conversely, a wide support or
#' \eqn{\nu \to 1} lets the errors absorb the moments and pulls \eqn{\hat P}
#' toward the prior.
#'
#' Standard errors of \eqn{\lambda} (\code{se_lambda}) use the analytic Hessian
#' of the dual (Golan 2008, eq. 4.7): the signal covariance
#' \eqn{\nu^{-1}\sum_k A_{ks}A_{ks'}P_{kj}(\delta_{jj'}-P_{kj'})} plus the noise
#' covariance \eqn{(1-\nu)^{-1}\sum_{i,t} z_{its} z_{its'}\mathrm{Var}_w(v)_{itj}}
#' (block-diagonal in the destination \eqn{j}). The noise term makes the Hessian
#' full rank -- \eqn{\lambda} is identified, so no reference normalization is
#' needed (unlike \code{\link{markov_ce}}). The reported SE is the robust
#' \strong{unit-clustered} sandwich \eqn{H^{-1}\hat V H^{-1}} (Golan 2008, Sec.
#' 3.3), which matches a unit bootstrap; the naive \eqn{H^{-1}} and the
#' non-clustered sandwich are conservative because covariate scores are correlated
#' within a unit.
#'
#' @param data A balanced long-format \code{data.frame}, one row per unit-period.
#' @param id,time Names of the unit-identifier and time columns.
#' @param states Character vector of state columns; each row over them is
#'   non-negative and sums to 1. \code{K = length(states)}.
#' @param covariates Character vector of covariate columns
#'   (\code{S = length(covariates)}); required. Time-varying covariates are
#'   supported.
#' @param v Numeric error support, symmetric about 0 with at least two points;
#'   default \code{c(-1, 0, 1)} (Golan's symmetric \eqn{[-1,1]} support).
#' @param nu Signal/noise weight in \eqn{(0,1)}; default \code{0.5}. Smaller
#'   \code{nu} penalizes the noise more strongly (driving \eqn{\hat\varepsilon
#'   \to 0}, toward the exact \code{markov_ce}); larger \code{nu} lets the noise
#'   absorb more of the moment discrepancy.
#' @param p0 Optional \code{K}-by-\code{K} prior transition matrix, strictly
#'   positive with rows summing to 1. Defaults to uniform.
#' @param w0 Optional error prior on the support (length \code{length(v)},
#'   positive, summing to 1). Defaults to uniform (mean-zero error prior).
#' @param control A named list merged over the defaults and passed to
#'   \code{\link[stats]{optim}} (BFGS): \code{maxit} (default 1000) and
#'   \code{reltol} (default 1e-12). \code{fnscale} is forced to -1 (the dual is
#'   maximized).
#'
#' @return An object of class \code{c("markov_gce", "infometrics")} with \code{p}/
#'   \code{p_hat} (\code{K x K} row-stochastic), \code{lambda}/\code{lambda_hat}
#'   (\code{S x K}), \code{epsilon} (fitted per-observation errors, Pn x K),
#'   \code{se} (multinomial SE of \code{P}), \code{se_lambda} (\code{S x K}
#'   analytic-Hessian unit-clustered sandwich SE of \eqn{\lambda}; all finite --
#'   the noise identifies \eqn{\lambda}, so no reference is needed),
#'   \code{vcov_lambda}, \code{states}, \code{n_transitions}, \code{nu}, \code{v},
#'   \code{H_p}/\code{H_signal} (signal/mean row entropy), \code{H_w} (noise
#'   entropy), \code{entropy}, \code{S} (normalized entropy),
#'   \code{fitted.values}, \code{residuals}, \code{moment_residual}
#'   (exact-moment gap, expected nonzero), \code{foc_residual} (noisy-moment gap,
#'   ~0 at the optimum), \code{objective}/\code{value}, \code{prior},
#'   \code{converged}, \code{convergence}, \code{method}, and \code{call}.
#'
#' @section Choosing nu and v: There is no universally correct default. The
#'   noise enters per observation and the literal support \eqn{v = (-1,0,1)} is
#'   permissive relative to share-valued moment contributions, so at a wide
#'   support or large \code{nu} the errors can absorb the signal and pull
#'   \eqn{\hat P} toward the prior. In controlled recovery experiments the
#'   estimate approaches the exact \code{\link{markov_ce}} solution monotonically
#'   as the support narrows (\eqn{v \to 0}) or as \code{nu} decreases toward 0,
#'   though an ultra-narrow support can make the over-identified exact problem
#'   numerically unbounded (non-convergence). Users should report sensitivity to
#'   \code{nu} and \code{v} for their data rather than rely on a single setting.
#'
#' @references Golan, A. (2008). \emph{Information and Entropy Econometrics -
#'   A Review and Synthesis}. Foundations and Trends in Econometrics, 2(1-2),
#'   1-145. Section 7.7.1, eqs. (7.25a), (7.28).
#' @seealso \code{\link{markov_ce}}, \code{\link{matrix_gce}}
#'
#' @examples
#' set.seed(1)
#' P_true <- matrix(c(.6,.3,.1, .2,.5,.3, .1,.3,.6), nrow = 3, byrow = TRUE)
#' N <- 120L; Tn <- 5L; K <- 3L
#' Y <- matrix(rgamma(N * K, 1), N, K); Y <- Y / rowSums(Y); inc <- rnorm(N)
#' id <- rep(seq_len(N), each = Tn); tm <- rep(seq_len(Tn), N)
#' sh <- matrix(0, N * Tn, K); income <- numeric(N * Tn); infl <- numeric(N * Tn)
#' for (t in seq_len(Tn)) {
#'   idx <- which(tm == t); inc <- 0.8 * inc + rnorm(N, 0, 0.5)
#'   sh[idx, ] <- Y; income[idx] <- inc; infl[idx] <- 0.02 * t
#'   Y <- (Y %*% P_true) * exp(matrix(rnorm(N * K, 0, 0.03), N, K)); Y <- Y / rowSums(Y)
#' }
#' panel <- data.frame(i = id, t = tm, y1 = sh[, 1], y2 = sh[, 2], y3 = sh[, 3],
#'                     income = income, infl = infl)
#' fit <- markov_gce(panel, id = "i", time = "t", states = c("y1","y2","y3"),
#'                   covariates = c("income", "infl"), v = c(-0.2, 0, 0.2),
#'                   nu = 0.2)
#' coef(fit)
#'
#' @importFrom stats optim sd pnorm
#' @export
markov_gce <- function(data, id, time, states, covariates,
                       v = NULL, nu = 0.5, p0 = NULL, w0 = NULL,
                       control = list()) {
  if (!is.data.frame(data)) stop("'data' must be a data.frame in long format.")
  if (length(states) < 2L) stop("'states' must name at least two state columns.")
  if (missing(covariates) || length(covariates) < 1L)
    stop("'covariates' is required for markov_gce(); use markov_ce() for the ",
         "covariate-free / exact case.")
  if (!is.numeric(nu) || length(nu) != 1L || nu <= 0 || nu >= 1)
    stop("'nu' must be a single number in (0, 1).")

  if (is.null(v)) v <- c(-1, 0, 1)
  v <- as.numeric(v)
  if (length(v) < 2L) stop("'v' must have at least two support points.")
  if (any(abs(v + rev(v)) > 1e-8)) stop("'v' must be symmetric about 0.")
  Mv <- length(v)
  if (is.null(w0)) w0 <- rep(1 / Mv, Mv)
  else { w0 <- as.numeric(w0)
         if (length(w0) != Mv || any(w0 <= 0)) stop("'w0' must be positive, length(v).")
         w0 <- w0 / sum(w0) }
  logw0 <- log(w0)

  a <- .markov_arrays(data, id, time, states, covariates)
  N <- a$N; Tn <- a$Tn; K <- a$K; lab <- a$states; S <- dim(a$Z)[3]

  flat <- function(arr, V) t(matrix(aperm(arr[, , , drop = FALSE], c(3, 1, 2)), nrow = V))
  XX <- flat(a$Y[, 1:(Tn - 1), , drop = FALSE], K)       # P x K  lagged states (P = N*(Tn-1))
  YY <- flat(a$Y[, 2:Tn,       , drop = FALSE], K)       # lead states
  Z2 <- flat(a$Z[, 1:(Tn - 1), , drop = FALSE], S)       # lagged covariates
  Z1 <- flat(a$Z[, 2:Tn,       , drop = FALSE], S)       # lead covariates
  Pn <- nrow(XX)

  ## standardize covariates for conditioning; report lambda on original scale
  cov_scale <- apply(rbind(Z1, Z2), 2L, function(z) { s <- stats::sd(z); if (s > 0) s else 1 })
  Z1s <- sweep(Z1, 2L, cov_scale, "/"); Z2s <- sweep(Z2, 2L, cov_scale, "/")

  ## signal cross-moments (same as markov_ce, on standardized covariates)
  A <- crossprod(XX,  Z2s)                               # K x S
  M <- crossprod(Z1s, YY)                                # S x K

  if (is.null(p0)) p0 <- matrix(1 / K, K, K) else {
    p0 <- as.matrix(p0)
    if (!all(dim(p0) == c(K, K))) stop("p0 must be K x K = ", K, " x ", K, ".")
    if (any(p0 <= 0)) stop("All prior masses 'p0' must be strictly positive.")
    rs <- rowSums(p0); if (any(abs(rs - 1) > 1e-8)) { warning("Rows of 'p0' renormalized."); p0 <- p0 / rs }
  }
  logp0 <- log(p0)

  ## ---- pieces: recover p (K x K) and the noise expectation eps (Pn x K) ----
  pieces <- function(lvec) {
    lambda <- matrix(lvec, S, K)
    ## signal
    eta <- (A %*% lambda) / nu                           # K x K
    Ls  <- logp0 + eta
    mk  <- apply(Ls, 1L, max); ELs <- exp(Ls - mk); rk <- rowSums(ELs)
    p   <- ELs / rk                                      # K x K
    logOmega <- mk + log(rk)
    ## noise (per observation r and destination j)
    theta <- (Z1s %*% lambda) / (1 - nu)                 # Pn x K
    logPsi <- matrix(0, Pn, K); eps <- matrix(0, Pn, K)
    for (j in seq_len(K)) {
      TM  <- outer(theta[, j], v) + matrix(logw0, Pn, Mv, byrow = TRUE)  # Pn x Mv
      mr  <- apply(TM, 1L, max); ETM <- exp(TM - mr); rr <- rowSums(ETM)
      logPsi[, j] <- mr + log(rr)
      what <- ETM / rr                                   # Pn x Mv
      eps[, j] <- what %*% v
    }
    list(lambda = lambda, p = p, logOmega = logOmega, eps = eps, logPsi = logPsi)
  }

  obj <- function(lvec) {
    pc <- pieces(lvec)
    sum(M * pc$lambda) - nu * sum(pc$logOmega) - (1 - nu) * sum(pc$logPsi)
  }
  ## gradient: dl/dlambda_sj = M_sj - sum_k A_ks p_kj - sum_r Z1s[r,s] eps_rj
  grad <- function(lvec) {
    pc <- pieces(lvec)
    as.vector(M - t(A) %*% pc$p - crossprod(Z1s, pc$eps))
  }

  ## maximise the dual (fnscale = -1, kept from the draft by user choice)
  con <- list(maxit = 1000L, reltol = 1e-12)
  con[names(control)] <- control
  con$fnscale <- -1
  est <- stats::optim(rep(0, S * K), fn = obj, gr = grad, method = "BFGS", control = con)
  if (est$convergence != 0L) warning("optim did not converge (code ", est$convergence, ").")

  pc <- pieces(est$par)
  p  <- pc$p; dimnames(p) <- list(lab, lab)
  lambda <- pc$lambda / cov_scale; rownames(lambda) <- covariates; colnames(lambda) <- lab

  n_from <- colSums(XX)
  se <- sqrt(p * (1 - p) / pmax(n_from, 1)); dimnames(se) <- list(lab, lab)
  fitted <- XX %*% p

  ## diagnostics
  ent_row <- apply(p,  1L, function(z) { z <- z[z > 0]; -sum(z * log(z)) })
  H0_row  <- apply(p0, 1L, function(z) { z <- z[z > 0]; -sum(z * log(z)) })
  Hp <- mean(ent_row); Hp0 <- mean(H0_row)
  ## noise entropy: recompute w to summarize H(w)
  Hw <- {
    theta <- (Z1s %*% pc$lambda) / (1 - nu); tot <- 0; cnt <- 0
    for (j in seq_len(K)) {
      TM <- outer(theta[, j], v) + matrix(logw0, Pn, Mv, byrow = TRUE)
      mr <- apply(TM, 1L, max); ETM <- exp(TM - mr); what <- ETM / rowSums(ETM)
      tot <- tot + sum(apply(what, 1L, function(w) { w <- w[w > 0]; -sum(w * log(w)) })); cnt <- cnt + Pn
    }
    tot / cnt
  }

  ## score inputs for margins() sandwich SEs (covariates always present here)
  se_data <- list(A = A, M = M, XX = XX, YY = YY, Z1s = Z1s, Z2s = Z2s,
                  cov_scale = cov_scale, logp0 = logp0, nu = nu, n_from = n_from,
                  v = v, logw0 = logw0, unit = rep(seq_len(N), Tn - 1L))

  ## analytic-Hessian (Golan eq. 4.7: signal/nu + block-diagonal noise/(1-nu))
  ## unit-clustered sandwich SEs of lambda -- full rank (the noise identifies
  ## lambda), so no reference normalization is needed.
  lse <- .markov_lambda_se_full(se_data, lambda)
  se_lambda <- lse$se; dimnames(se_lambda) <- list(covariates, lab)

  structure(
    list(p = p, p_hat = p, lambda = lambda, lambda_hat = lambda,
         epsilon = pc$eps, se = se, se_lambda = se_lambda,
         vcov_lambda = lse$vcov, states = lab, se_data = se_data,
         n_transitions = Pn, n_from = n_from, nu = nu, v = v, w0 = w0,
         H_p = Hp, H_signal = Hp, H_w = Hw, entropy = Hp,
         S = if (Hp0 > 0) Hp / Hp0 else NA_real_,
         fitted.values = fitted, residuals = YY - fitted,
         moment_residual = max(abs(M - t(A) %*% p)),                  # exact-moment gap (nonzero)
         foc_residual = max(abs(M - t(A) %*% p - crossprod(Z1s, pc$eps))),  # noisy FOC (~0)
         objective = est$value, value = est$value, prior = p0,
         converged = (est$convergence == 0L), convergence = est$convergence,
         method = "dual", call = match.call()),
    class = c("markov_gce", "infometrics")
  )
}

## ---- S3 methods ------------------------------------------------------------
#' @rdname markov_gce
#' @param object,x A \code{markov_gce} object.
#' @param ... Unused.
#' @export
coef.markov_gce <- function(object, ...) object$p
#' @rdname markov_gce
#' @export
fano_bounds.markov_gce <- function(object, ...)
  .fano_row_bounds(object$p_hat, ncol(object$p_hat))
#' @rdname markov_gce
#' @param average Logical; if \code{TRUE} (default) return the origin-weighted
#'   average marginal effect (\code{S x K}); if \code{FALSE} the full
#'   \code{K x K x S} array of per-origin effects.
#' @param se Logical (default \code{FALSE}); if \code{TRUE} (with
#'   \code{average = TRUE}) attach standard errors.
#' @param se_method \code{"sandwich"} (default) or \code{"bootstrap"}. The naive
#'   Hessian-inverse delta overstates ~4-5x here; the robust sandwich
#'   \eqn{(-H)^{-1}\hat V(-H)^{-1}} matches the \code{"bootstrap"} (which resamples
#'   units and re-solves the dual).
#' @param B Integer bootstrap resamples when \code{se_method = "bootstrap"}
#'   (default 500).
#' @export
margins.markov_gce <- function(object, average = TRUE, se = FALSE,
                               se_method = c("sandwich", "bootstrap"), B = 500L, ...)
  .margins_markov(object, average, se, match.arg(se_method), B)
#' @rdname markov_gce
#' @export
fitted.markov_gce <- function(object, ...) object$fitted.values
#' @rdname markov_gce
#' @export
residuals.markov_gce <- function(object, ...) object$residuals
#' @rdname markov_gce
#' @param digits Significant digits to print.
#' @export
print.markov_gce <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Markov matrix via Generalized Cross-Entropy (", length(x$states),
      " states, ", x$n_transitions, " transitions, nu = ", x$nu, ")\n", sep = "")
  cat("\nP (rows sum to 1, P[k, j] = k -> j):\n"); print.default(round(x$p, digits))
  cat("\nH(p) =", format(x$H_p, digits = digits), " H(w) =", format(x$H_w, digits = digits),
      " normalized S =", format(x$S, digits = digits), "\n")
  cat("noisy-moment FOC residual =", format(x$foc_residual, digits = digits),
      " (exact-moment gap =", format(x$moment_residual, digits = digits), ")\n")
  cat("convergence:", x$convergence,
      if (x$convergence == 0L) "(converged)" else "(see ?optim)", "\n")
  invisible(x)
}
#' @rdname markov_gce
#' @export
summary.markov_gce <- function(object, ...) {
  rows <- apply(object$p, 1L, function(z) { z <- z[z > 0]; -sum(z * log(z)) })
  structure(list(call = object$call, p = object$p, se = object$se, lambda = object$lambda,
                 se_lambda = object$se_lambda,
                 states = object$states, n_transitions = object$n_transitions,
                 nu = object$nu, v = object$v, H_p = object$H_p, H_w = object$H_w,
                 row_S = rows / log(length(object$states)), entropy = object$entropy,
                 S = object$S, foc_residual = object$foc_residual,
                 moment_residual = object$moment_residual, convergence = object$convergence),
            class = "summary.markov_gce")
}
#' @export
print.summary.markov_gce <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("GCE transition matrix P (rows sum to 1), nu = ", x$nu,
      ", support v = (", paste(x$v, collapse = ", "), "):\n", sep = "")
  print.default(round(x$p, digits))
  cat("\nCovariate effects lambda (S x K, original scale):\n"); print.default(round(x$lambda, digits))
  if (!is.null(x$se_lambda) && any(is.finite(x$se_lambda))) {
    cat("\nSE(lambda) [analytic-Hessian unit-clustered sandwich, Golan eq. 4.7 / sec 3.3]:\n")
    print.default(round(x$se_lambda, digits))
  }
  cat("\nPer-row normalized entropy S:", paste(format(x$row_S, digits = digits), collapse = "  "), "\n")
  fb <- attr(.fano_row_bounds(x$p, ncol(x$p)), "overall")
  cat("Fano (sec 7.5): mean next-state error pe =", format(fb["mean_pe"], digits = digits),
      ">= bound", format(fb["mean_pe_lower"], digits = digits), "\n")
  cat("H(p) =", format(x$H_p, digits = digits), " H(w) =", format(x$H_w, digits = digits),
      " normalized S(P) =", format(x$S, digits = digits), "\n")
  cat("noisy-moment FOC residual =", format(x$foc_residual, digits = digits),
      " exact-moment gap =", format(x$moment_residual, digits = digits), "\n")
  cat("convergence:", x$convergence,
      if (x$convergence == 0L) "(converged)" else "(see ?optim)", "\n")
  invisible(x)
}
