# markov_ce.R
# Cross-Entropy estimation of a first-order, stationary Markov transition /
# balancing matrix from a balanced long panel (Golan 2008, Section 7.7.1,
# eqs. 7.25-7.27). States are simplex-valued (one-hot indicators OR
# compositional shares). Optional covariates z enter lead/lagged like the
# states: per unit and period, Z2 = z(t-1) (with the lagged shares),
# Z1 = z(t) (with the lead shares).
#   A = sum_{i,t} XX' Z2   (K x S);  M = sum_{i,t} Z1' YY   (S x K)
#   lambda is S x K;  dual  l(lambda) = sum(M*lambda) - sum_k log Omega_k.
# Default (no covariates) uses the lagged-state indicator (Case A).
#
# Deliberate standalone estimator (panel/time-indexed cousin of matrix_ce).
# By user instruction it keeps its own argument names (data, id, time, states,
# covariates, p0) rather than me()'s. The dual is kept in the draft's
# maximise form (fnscale = -1) by user choice -- the one estimator in the
# package that maximises rather than minimises.
#
# References:
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145.
#     Section 7.7.1.

## ---- internal: reshape a balanced long panel to N x Tn x V arrays ----------

#' @keywords internal
.markov_arrays <- function(data, id, time, states, covariates) {
  for (cl in c(id, time, states, covariates))
    if (!cl %in% names(data)) stop("Column '", cl, "' not found in 'data'.")

  units <- unique(data[[id]]); periods <- sort(unique(data[[time]]))
  N <- length(units); Tn <- length(periods); K <- length(states)
  if (N * Tn != nrow(data))
    stop("Panel is not balanced: N*T = ", N * Tn, " but nrow(data) = ", nrow(data),
         ". markov_ce() requires a balanced panel (every id at every time).")
  if (Tn < 2L) stop("Need at least two periods to form a transition.")

  d   <- data[order(match(data[[id]], units), match(data[[time]], periods)), , drop = FALSE]
  tab <- table(factor(d[[id]], levels = units), factor(d[[time]], levels = periods))
  if (any(tab != 1L))
    stop("Panel is not balanced (missing or duplicated periods), e.g. unit '",
         units[which(rowSums(tab != 1L) > 0)[1]], "'.")

  Ymat <- as.matrix(d[, states, drop = FALSE])
  if (any(!is.finite(Ymat)) || any(Ymat < 0))
    stop("State columns must be non-negative and non-missing.")
  if (any(abs(rowSums(Ymat) - 1) > 1e-6))
    stop("Each row's state values (", paste(states, collapse = ", "),
         ") must sum to 1 (one-hot indicators or compositional shares).")
  Yarr <- aperm(array(Ymat, dim = c(Tn, N, K)), c(2, 1, 3))   # N x Tn x K
  dimnames(Yarr) <- list(NULL, as.character(periods), states)

  Zarr <- NULL
  if (!is.null(covariates)) {
    Zmat <- as.matrix(d[, covariates, drop = FALSE])
    if (any(!is.finite(Zmat))) stop("Covariate columns must be finite and non-missing.")
    S <- length(covariates)
    Zarr <- aperm(array(Zmat, dim = c(Tn, N, S)), c(2, 1, 3))  # N x Tn x S
    dimnames(Zarr) <- list(NULL, as.character(periods), covariates)
  }
  list(Y = Yarr, Z = Zarr, N = N, Tn = Tn, K = K, states = states)
}

## ---- shared internals for fano_bounds() and margins() (markov_ce + markov_gce)

#' @keywords internal
.fd_jac <- function(fn, x, h = 1e-5) {
  n <- length(x); f0 <- fn(x); J <- matrix(0, length(f0), n)
  for (i in seq_len(n)) { e <- numeric(n); e[i] <- h; J[, i] <- (fn(x + e) - fn(x - e)) / (2 * h) }
  J
}

# P (K x K) and, for the noisy sibling, the per-transition error eps (Pn x K),
# as a function of the standardized multiplier vector; reused by the sandwich.
#' @keywords internal
.markov_P_eps <- function(lvec, sd) {
  Kk <- nrow(sd$A); S <- ncol(sd$A); nu <- sd$nu
  lam <- matrix(lvec, S, Kk)
  Ls  <- sd$logp0 + (sd$A %*% lam) / nu
  mk  <- apply(Ls, 1L, max); EL <- exp(Ls - mk); rk <- rowSums(EL); P <- EL / rk
  logOmega <- mk + log(rk)                               # K
  eps <- NULL; logPsi <- NULL
  if (!is.null(sd$v)) {                                   # gce: noise term present
    Mv <- length(sd$v); Pn <- nrow(sd$Z1s); theta <- (sd$Z1s %*% lam) / (1 - nu)
    eps <- matrix(0, Pn, Kk); logPsi <- matrix(0, Pn, Kk)
    for (j in seq_len(Kk)) {
      TM <- outer(theta[, j], sd$v) + matrix(sd$logw0, Pn, Mv, byrow = TRUE)
      mr <- apply(TM, 1L, max); ETM <- exp(TM - mr); rr <- rowSums(ETM)
      logPsi[, j] <- mr + log(rr); eps[, j] <- (ETM / rr) %*% sd$v
    }
  }
  list(P = P, eps = eps, logOmega = logOmega, logPsi = logPsi)
}

# Re-solve the (maximised) dual for the standardized multipliers given a design
# list (used by the unit bootstrap). Unified over ce (nu = 1, no noise) and gce.
#' @keywords internal
.markov_refit <- function(sd, control = list()) {
  Kk <- nrow(sd$A); S <- ncol(sd$A)
  fn <- function(lvec) {
    pc <- .markov_P_eps(lvec, sd)
    val <- sum(sd$M * matrix(lvec, S, Kk)) - sd$nu * sum(pc$logOmega)
    if (!is.null(pc$logPsi)) val <- val - (1 - sd$nu) * sum(pc$logPsi)
    val
  }
  gr <- function(lvec) {
    pc <- .markov_P_eps(lvec, sd)
    g  <- sd$M - crossprod(sd$A, pc$P)
    if (!is.null(pc$eps)) g <- g - crossprod(sd$Z1s, pc$eps)
    as.vector(g)
  }
  con <- list(maxit = 1000L, reltol = 1e-12); con[names(control)] <- control
  con$fnscale <- -1
  est <- stats::optim(rep(0, S * Kk), fn, gr, method = "BFGS", control = con)
  list(par = est$par, converged = est$convergence == 0L)
}

# Unit (block) bootstrap SE of the AME: resample units, rebuild the transition
# design from se_data, re-solve, recompute the AME. Standardization cancels, so
# this equals a refit-from-data bootstrap.
#' @keywords internal
.markov_boot_se <- function(sd, S, Kk, nu, B, dn) {
  units <- unique(sd$unit); rows_of <- split(seq_along(sd$unit), sd$unit)
  AMEb <- array(NA_real_, c(S, Kk, B))
  for (b in seq_len(B)) {
    rows <- unlist(rows_of[as.character(sample(units, length(units), replace = TRUE))],
                   use.names = FALSE)
    sdb <- sd
    sdb$XX <- sd$XX[rows, , drop = FALSE]; sdb$YY <- sd$YY[rows, , drop = FALSE]
    sdb$Z1s <- sd$Z1s[rows, , drop = FALSE]; sdb$Z2s <- sd$Z2s[rows, , drop = FALSE]
    sdb$A <- crossprod(sdb$XX, sdb$Z2s); sdb$M <- crossprod(sdb$Z1s, sdb$YY)
    sdb$n_from <- colSums(sdb$XX)
    rf <- tryCatch(.markov_refit(sdb), error = function(e) NULL)
    if (is.null(rf) || !rf$converged) next
    pc <- .markov_P_eps(rf$par, sdb)
    AMEb[, , b] <- .markov_margins(pc$P, matrix(rf$par, S, Kk) / sd$cov_scale,
                                   nu, sdb$n_from, average = TRUE)
  }
  n_ok <- sum(!is.na(AMEb[1, 1, ]))
  if (n_ok < 2L) stop("bootstrap produced < 2 usable resamples.")
  SE <- matrix(apply(AMEb, c(1L, 2L), function(z) stats::sd(z, na.rm = TRUE)), S, Kk,
               dimnames = dn)
  attr(SE, "n_boot") <- n_ok
  SE
}

# Analytic Hessian of the CE log-partition (Golan 2008, eq. 4.7): the covariance
# of the moment functions under the fitted P.  H_{(s,j),(s',j')} =
# sum_k A_ks A_ks' P_kj(delta_jj' - P_kj')  (vec order: s fast within j).  This is
# the -d(grad)/d(lambda) that .markov_sandwich otherwise finite-differences.
#' @keywords internal
.markov_hessian <- function(A, P) {
  Kk <- nrow(A); S <- ncol(A); H <- matrix(0, S * Kk, S * Kk)
  for (j in seq_len(Kk)) for (jp in seq_len(Kk)) {
    w   <- P[, j] * ((j == jp) - P[, jp])                 # length K (over origins)
    blk <- crossprod(A, A * w)                            # S x S: sum_k w_k a_k a_k'
    H[((j - 1) * S + 1):(j * S), ((jp - 1) * S + 1):(jp * S)] <- blk
  }
  H
}

# Full analytic Hessian (information matrix) of the dual, unified over the exact
# CE (markov_ce: nu = 1, no noise) and the noisy GCE (markov_gce). Signal part
# = .markov_hessian(A, P) / nu; the noise part (eq. 4.7 for the noise partition
# Psi_itj, which depends only on column j of lambda) is block-diagonal in the
# destination j: (1/(1-nu)) sum_{i,t} z1_its z1_its' Var_w(v)_itj.
#' @keywords internal
.markov_info <- function(sd, lhat) {
  Kk <- nrow(sd$A); S <- ncol(sd$A)
  H  <- .markov_hessian(sd$A, .markov_P_eps(lhat, sd)$P) / sd$nu
  if (!is.null(sd$v)) {                                    # gce: add the noise term
    Pn <- nrow(sd$Z1s); Mv <- length(sd$v)
    theta <- (sd$Z1s %*% matrix(lhat, S, Kk)) / (1 - sd$nu)
    for (j in seq_len(Kk)) {
      TM <- outer(theta[, j], sd$v) + matrix(sd$logw0, Pn, Mv, byrow = TRUE)
      w  <- exp(TM - apply(TM, 1L, max)); w <- w / rowSums(w)
      e  <- w %*% sd$v; Varw <- as.vector(w %*% (sd$v^2) - e^2)
      idx <- ((j - 1) * S + 1):(j * S)
      H[idx, idx] <- H[idx, idx] + crossprod(sd$Z1s, sd$Z1s * Varw) / (1 - sd$nu)
    }
  }
  H
}

# Reference-normalized standard errors of lambda for the exact CE (markov_ce).
# Golan eq. 4.7 gives the Hessian H, but H^{-1} alone OVERSTATES the SE by ~sqrt
# (n_from) here because the moments A, M are cross-moment sums (H = n_from * Vhat);
# the correct SE is the sandwich H^{-1} Vhat H^{-1} (Golan Sec. 3.3), which matches
# a unit bootstrap.  lambda is identified only up to a per-conditioning additive
# shift, so SEs are relative to a reference destination `ref` (that column is NA);
# a boundary fit (P_kj in {0,1}) makes H singular -> all-NA (tryCatch). Returns SEs
# on the original covariate scale (rescaled by 1 / cov_scale when standardized).
#' @keywords internal
.markov_lambda_se <- function(A, P, XX, YY, Z1s, Z2s, cov_scale, ref) {
  Kk <- nrow(A); S <- ncol(A)
  H  <- .markov_hessian(A, P)
  fitted <- XX %*% P
  G <- matrix(0, nrow(XX), S * Kk)                        # per-transition scores
  for (j in seq_len(Kk)) G[, ((j - 1) * S + 1):(j * S)] <- Z1s * YY[, j] - Z2s * fitted[, j]
  Vhat <- crossprod(G)
  free <- which(rep(seq_len(Kk), each = S) != ref)        # drop reference destination
  se <- matrix(NA_real_, S, Kk); vc_out <- NULL
  vc <- tryCatch({
    Hi <- solve(H[free, free, drop = FALSE])
    Hi %*% Vhat[free, free, drop = FALSE] %*% Hi
  }, error = function(e) NULL)
  if (!is.null(vc) && all(is.finite(vc))) {
    if (!is.null(cov_scale)) {                            # standardized -> original scale
      D  <- 1 / cov_scale[rep(seq_len(S), Kk)[free]]
      vc <- outer(D, D) * vc
    }
    vc_out <- vc
    se[, -ref] <- matrix(sqrt(pmax(diag(vc), 0)), S, Kk - 1L)
  }
  list(se = se, vcov = vc_out, ref = ref)
}

# Full (no reference) UNIT-CLUSTERED analytic-Hessian sandwich SE of lambda for
# the noisy GCE (markov_gce): its Hessian is full rank (the noise breaks the
# softmax shift-redundancy, so lambda is identified). Clustering by unit accounts
# for the within-unit dependence (the non-clustered sandwich overstates ~1.26x;
# the clustered one matches a unit bootstrap). Returns an S x K SE (no NA column)
# on the original covariate scale; all-NA on a singular / boundary fit.
#' @keywords internal
.markov_lambda_se_full <- function(sd, lambda) {
  Kk <- nrow(sd$A); S <- ncol(sd$A)
  lhat <- as.vector(lambda * sd$cov_scale)
  H  <- .markov_info(sd, lhat)
  pc <- .markov_P_eps(lhat, sd); fitted <- sd$XX %*% pc$P
  G  <- matrix(0, nrow(sd$XX), S * Kk)
  for (j in seq_len(Kk)) {
    yj <- sd$YY[, j]; if (!is.null(pc$eps)) yj <- yj - pc$eps[, j]
    G[, ((j - 1) * S + 1):(j * S)] <- sd$Z1s * yj - sd$Z2s * fitted[, j]
  }
  Gc <- rowsum(G, sd$unit)                                # sum scores within each unit
  se <- matrix(NA_real_, S, Kk); vc_out <- NULL
  vc <- tryCatch({ Hi <- solve(H); Hi %*% crossprod(Gc) %*% Hi }, error = function(e) NULL)
  if (!is.null(vc) && all(is.finite(vc))) {
    D  <- 1 / sd$cov_scale[rep(seq_len(S), Kk)]           # standardized -> original scale
    vc <- outer(D, D) * vc
    vc_out <- vc
    se[] <- sqrt(pmax(diag(vc), 0))
  }
  list(se = se, vcov = vc_out)
}

# Robust sandwich (-H)^{-1} Vhat (-H)^{-1} for the standardized multipliers, with
# Vhat = sum_r g_r g_r' the per-transition dual scores
# g_r,sj = z1_rs (y_rj - eps_rj) - z2_rs (x_r'P)_j.  The naive Hessian-inverse
# overstates ~4-5x for the noisy Markov model; the sandwich matches a unit bootstrap.
#' @keywords internal
.markov_sandwich <- function(object) {
  sd <- object$se_data; Kk <- nrow(sd$A); S <- ncol(sd$A)
  lhat <- as.vector(object$lambda * sd$cov_scale)         # internal (standardized) lambda
  Ainv <- solve(.markov_info(sd, lhat))                   # (-H)^{-1}, analytic bread (eq. 4.7)
  pc   <- .markov_P_eps(lhat, sd); fitted <- sd$XX %*% pc$P
  Pn   <- nrow(sd$XX); G <- matrix(0, Pn, S * Kk)
  for (r in seq_len(Pn)) {
    yr <- sd$YY[r, ]; if (!is.null(pc$eps)) yr <- yr - pc$eps[r, ]
    G[r, ] <- as.vector(outer(sd$Z1s[r, ], yr) - outer(sd$Z2s[r, ], fitted[r, ]))
  }
  Ainv %*% crossprod(G) %*% Ainv
}

# Marginal effects of covariates on the transition probabilities:
# ME[k,j,s] = (P_kj/nu)(lambda_sj - sum_l P_kl lambda_sl). average = TRUE returns
# the origin-frequency-weighted AME (S x K); FALSE the full K x K x S array.
#' @keywords internal
.markov_margins <- function(P, lambda, nu, n_from, average) {
  Kk <- ncol(P); S <- nrow(lambda)
  ME <- array(0, c(Kk, Kk, S),
              dimnames = list(rownames(P), colnames(P), rownames(lambda)))
  for (k in seq_len(Kk)) {
    Pk <- P[k, ]; Jk <- diag(Pk) - outer(Pk, Pk)
    for (s in seq_len(S)) ME[k, , s] <- (Jk %*% lambda[s, ]) / nu
  }
  if (!average) return(ME)
  pw  <- n_from / sum(n_from)
  AME <- matrix(0, S, Kk, dimnames = list(rownames(lambda), colnames(P)))
  for (s in seq_len(S)) AME[s, ] <- colSums(ME[, , s] * pw)   # sum_k pw_k ME[k,,s]
  AME
}

# Shared margins() body: point estimates, or a margins_gce object with SEs
# (robust sandwich, default; or unit bootstrap).
#' @keywords internal
.margins_markov <- function(object, average, se, se_method = "sandwich", B = 500L) {
  if (is.null(object$se_data))
    stop("marginal effects require covariates; refit with 'covariates = ' ",
         "(the covariate-free transition matrix has nothing to differentiate).")
  sd <- object$se_data; nu <- sd$nu
  ME <- .markov_margins(object$p_hat, object$lambda, nu, sd$n_from, average)
  if (!se) return(ME)
  if (!average)
    stop("standard errors are implemented only for the averaged effects ",
         "(average = TRUE).")
  S <- nrow(object$lambda); Kk <- ncol(object$p_hat)
  if (se_method == "sandwich") {
    Vlam <- .markov_sandwich(object)                      # sandwich vcov of std lambda
    lhat <- as.vector(object$lambda * sd$cov_scale)
    me_ame <- function(lvec)
      as.vector(.markov_margins(.markov_P_eps(lvec, sd)$P, matrix(lvec, S, Kk) / sd$cov_scale,
                                nu, sd$n_from, average = TRUE))
    Jac <- .fd_jac(me_ame, lhat)
    SE  <- matrix(sqrt(pmax(diag(Jac %*% Vlam %*% t(Jac)), 0)), S, Kk,
                  dimnames = dimnames(ME))
  } else {                                                # se_method == "bootstrap"
    SE <- .markov_boot_se(sd, S, Kk, nu, B, dimnames(ME))
  }
  z <- ME / SE
  structure(list(estimate = ME, se = SE, z = z, p.value = 2 * stats::pnorm(-abs(z)),
                 se_method = se_method, n_boot = attr(SE, "n_boot")),
            class = "margins_gce")
}

#' Cross-Entropy Estimation of a Markov Transition/Balancing Matrix (Panel)
#'
#' Estimates a first-order, stationary row-stochastic matrix \eqn{P}
#' (\eqn{P_{kj} = k \to j}) linking simplex-valued state vectors over time,
#' from a balanced long panel, by Cross-Entropy (Golan 2008, Section 7.7.1,
#' eqs. 7.25-7.27). States may be one-hot membership indicators or
#' compositional shares (a panel/time-indexed extension of matrix balancing).
#' Data are supplied in long form (one row per unit-period) and reshaped
#' internally; transitions \eqn{t-1 \to t} are formed within each unit.
#'
#' @details
#' Optional covariates \code{z} enter lead/lagged like the states: per
#' transition the lagged states \eqn{X} pair with lagged covariates \eqn{Z_2}
#' and the lead states \eqn{Y} with lead covariates \eqn{Z_1}, forming the
#' cross moments \eqn{A = \sum X' Z_2} (\code{K x S}) and \eqn{M = \sum Z_1' Y}
#' (\code{S x K}), \eqn{S = }number of covariates. The concentrated dual
#' \deqn{\ell(\lambda) = \sum_{s,j} M_{sj}\lambda_{sj} - \sum_k \log\Omega_k(\lambda),
#'   \quad \hat P_{kj} = p0_{kj}\exp(\textstyle\sum_s A_{ks}\lambda_{sj})/\Omega_k}
#' is maximized over the \code{S x K} multipliers. With no covariates the
#' conditioning is the lagged-state indicator (\eqn{S = K}, transition counts).
#'
#' When the conditioning is full rank (\code{A = X'X}, the common indicator
#' case), the moment condition \eqn{A'P = M} exactly determines \eqn{P} as the
#' empirical/regression transition matrix, so the prior is irrelevant for
#' identified origin states. Only \strong{under-identified} rows
#' (\eqn{\mathrm{rank}(A) < K}, e.g. a state never observed as an origin) fall
#' back to the prior (maximum entropy when \code{p0} is uniform); see the
#' rank warning and the normalized entropy \code{S} in \code{summary()}.
#'
#' Standard errors of \eqn{\lambda} (\code{se_lambda}) use the analytic Hessian
#' of the CE log-partition (Golan 2008, eq. 4.7),
#' \eqn{H_{(s,j),(s',j')} = \sum_k A_{ks} A_{ks'} P_{kj}(\delta_{jj'} - P_{kj'})}
#' -- the covariance of the moment functions under \eqn{\hat P}. Note that
#' \eqn{H^{-1}} \emph{alone} overstates the SE (~\eqn{\sqrt{n}}), because the
#' moments \eqn{A, M} are cross-moment sums (\eqn{H = n\,\hat V}); the reported SE
#' is the robust sandwich \eqn{H^{-1}\hat V H^{-1}} (Golan 2008, Sec. 3.3), which
#' matches a unit bootstrap. Because \eqn{\lambda} is identified only up to a
#' per-conditioning additive shift, \code{se_lambda} is relative to a reference
#' destination (\code{lambda_ref = 1}, that column \code{NA}); a (near-)boundary
#' fit (\eqn{P_{kj}\in\{0,1\}}) yields \code{NA} (an infinite log-odds SE).
#'
#' @param data A balanced long-format \code{data.frame}, one row per unit-period.
#' @param id,time Names of the unit-identifier and time columns.
#' @param states Character vector of state columns; each row over them is
#'   non-negative and sums to 1. \code{K = length(states)}.
#' @param covariates Optional character vector of covariate columns
#'   (\code{S = length(covariates)}); default uses the lagged-state indicator.
#' @param p0 Optional \code{K}-by-\code{K} prior transition matrix, strictly
#'   positive with rows summing to 1. Defaults to uniform.
#' @param control A named list merged over the defaults and passed to
#'   \code{\link[stats]{optim}} (BFGS): \code{maxit} (default 1000) and
#'   \code{reltol} (default 1e-12). \code{fnscale} is forced to -1 (the dual is
#'   maximized).
#'
#' @return An object of class \code{c("markov_ce", "infometrics")} with \code{p}/
#'   \code{p_hat} (\code{K x K} row-stochastic), \code{lambda}/\code{lambda_hat}
#'   (\code{S x K}), \code{se} (multinomial SE of \code{P}), \code{se_lambda}
#'   (\code{S x K} analytic-Hessian sandwich SE of \eqn{\lambda}, reference
#'   destination column \code{NA}), \code{lambda_ref}, \code{vcov_lambda},
#'   \code{n_from}, \code{states},
#'   \code{n_transitions}, \code{fitted.values}, \code{residuals},
#'   \code{entropy}/\code{H_signal} (mean row entropy), \code{S} (prior-relative
#'   normalized entropy), \code{moment_residual}, \code{objective}/\code{value},
#'   \code{prior}, \code{state_type}, \code{converged}, \code{convergence},
#'   \code{method}, and \code{call}.
#'
#' @references Golan, A. (2008). \emph{Information and Entropy Econometrics -
#'   A Review and Synthesis}. Foundations and Trends in Econometrics, 2(1-2),
#'   1-145. Section 7.7.1.
#' @seealso \code{\link{matrix_ce}}, \code{\link{matrix_gce}}
#'
#' @examples
#' set.seed(1)
#' P_true <- matrix(c(.6,.3,.1, .2,.5,.3, .1,.3,.6), nrow = 3, byrow = TRUE)
#' N <- 150L; Tn <- 5L; K <- 3L
#' s  <- sample(K, N, replace = TRUE)
#' id <- rep(seq_len(N), each = Tn); tm <- rep(seq_len(Tn), N)
#' state <- integer(N * Tn)
#' for (t in seq_len(Tn)) {
#'   state[tm == t] <- s
#'   s <- vapply(s, function(k) sample(K, 1L, prob = P_true[k, ]), integer(1))
#' }
#' oh <- diag(K)[state, ]
#' panel <- data.frame(i = id, t = tm, y1 = oh[, 1], y2 = oh[, 2], y3 = oh[, 3])
#' fit <- markov_ce(panel, id = "i", time = "t", states = c("y1", "y2", "y3"))
#' coef(fit)
#'
#' @importFrom stats optim sd pnorm
#' @export
markov_ce <- function(data, id, time, states, covariates = NULL,
                      p0 = NULL, control = list()) {
  if (!is.data.frame(data)) stop("'data' must be a data.frame in long format.")
  if (length(states) < 2L) stop("'states' must name at least two state columns.")

  a <- .markov_arrays(data, id, time, states, covariates)
  N <- a$N; Tn <- a$Tn; K <- a$K; lab <- a$states

  ## within-unit lead/lag, stacked over transitions
  flat <- function(arr, V) t(matrix(aperm(arr[, , , drop = FALSE], c(3, 1, 2)), nrow = V))
  XX <- flat(a$Y[, 1:(Tn - 1), , drop = FALSE], K)       # (N*(Tn-1)) x K  lagged states
  YY <- flat(a$Y[, 2:Tn,       , drop = FALSE], K)       # lead states

  state_type <- if (all(abs(XX - round(XX)) < 1e-8)) "indicator" else "shares"

  if (is.null(covariates)) {
    Z2 <- XX; Z1 <- NULL                                  # Case A: lagged-state indicator
    A  <- crossprod(XX)                                   # K x K
    M  <- crossprod(XX, YY)                               # K x K
    S  <- K; cov_scale <- NULL
  } else {
    Sv <- dim(a$Z)[3]
    Z2 <- flat(a$Z[, 1:(Tn - 1), , drop = FALSE], Sv)     # lagged covariates
    Z1 <- flat(a$Z[, 2:Tn,       , drop = FALSE], Sv)     # lead covariates
    cov_scale <- apply(rbind(Z1, Z2), 2L, function(z) { s <- stats::sd(z); if (s > 0) s else 1 })
    Z1 <- sweep(Z1, 2L, cov_scale, "/"); Z2 <- sweep(Z2, 2L, cov_scale, "/")
    A  <- crossprod(XX, Z2)                               # K x S
    M  <- crossprod(Z1, YY)                               # S x K
    S  <- Sv
  }
  n_from <- colSums(XX)

  rankA <- qr(A)$rank
  if (rankA < K)
    warning("rank(A) = ", rankA, " < K = ", K, ": the conditioning does not ",
            "separate all origin states; the transition matrix is under-",
            "identified and unidentified rows fall back to the prior / maximum ",
            "entropy. See the normalized entropy S in summary().")

  if (is.null(p0)) {
    p0 <- matrix(1 / K, K, K)
  } else {
    p0 <- as.matrix(p0)
    if (!all(dim(p0) == c(K, K))) stop("p0 must be K x K = ", K, " x ", K, ".")
    if (any(p0 <= 0)) stop("All prior masses 'p0' must be strictly positive.")
    rs <- rowSums(p0)
    if (any(abs(rs - 1) > 1e-8)) { warning("Rows of 'p0' renormalized to sum to 1."); p0 <- p0 / rs }
  }
  logp0 <- log(p0)

  pieces <- function(lvec) {
    lambda <- matrix(lvec, S, K)
    L  <- logp0 + A %*% lambda
    mk <- apply(L, 1L, max); EL <- exp(L - mk); rk <- rowSums(EL)
    list(logOmega = mk + log(rk), p = EL / rk, lambda = lambda)
  }
  obj  <- function(lvec) { pc <- pieces(lvec); sum(M * pc$lambda) - sum(pc$logOmega) }
  grad <- function(lvec) { pc <- pieces(lvec); as.vector(M - t(A) %*% pc$p) }

  ## maximise the dual (fnscale = -1, kept from the draft by user choice)
  con <- list(maxit = 1000L, reltol = 1e-12)
  con[names(control)] <- control
  con$fnscale <- -1
  est <- stats::optim(rep(0, S * K), fn = obj, gr = grad, method = "BFGS", control = con)
  if (est$convergence != 0L) warning("optim did not converge (code ", est$convergence, ").")

  pc <- pieces(est$par)
  p  <- pc$p; dimnames(p) <- list(lab, lab)
  lambda <- pc$lambda
  if (!is.null(cov_scale)) {
    rownames(lambda) <- covariates
    lambda <- lambda / cov_scale          # report on the original covariate scale
  }

  se <- sqrt(p * (1 - p) / pmax(n_from, 1)); dimnames(se) <- list(lab, lab)
  fitted <- XX %*% p

  ## analytic-Hessian sandwich SEs of lambda (Golan eq. 4.7 bread), reference-
  ## normalized to destination `lambda_ref`; NA where the fit is (near-)boundary.
  lambda_ref <- 1L
  Z1_se <- if (is.null(covariates)) XX else Z1
  Z2_se <- if (is.null(covariates)) XX else Z2
  lse   <- .markov_lambda_se(A, p, XX, YY, Z1_se, Z2_se, cov_scale, lambda_ref)
  se_lambda <- lse$se
  dimnames(se_lambda) <- list(if (is.null(covariates)) lab else covariates, lab)

  ent_row <- apply(p,  1L, function(z) { z <- z[z > 0]; -sum(z * log(z)) })
  H0_row  <- apply(p0, 1L, function(z) { z <- z[z > 0]; -sum(z * log(z)) })
  Hp <- mean(ent_row); Hp0 <- mean(H0_row)

  ## score inputs for margins() sandwich SEs (only when covariates are present)
  se_data <- if (!is.null(cov_scale))
    list(A = A, M = M, XX = XX, YY = YY, Z1s = Z1, Z2s = Z2, cov_scale = cov_scale,
         logp0 = logp0, nu = 1, n_from = n_from, v = NULL, logw0 = NULL,
         unit = rep(seq_len(N), Tn - 1L)) else NULL

  structure(
    list(p = p, p_hat = p, lambda = lambda, lambda_hat = lambda,
         se = se, se_lambda = se_lambda, lambda_ref = lambda_ref,
         vcov_lambda = lse$vcov, n_from = n_from, states = lab, se_data = se_data,
         n_transitions = nrow(XX), state_type = state_type,
         fitted.values = fitted, residuals = YY - fitted,
         entropy = Hp, H_signal = Hp,
         S = if (Hp0 > 0) Hp / Hp0 else NA_real_,
         moment_residual = max(abs(M - t(A) %*% p)),
         objective = est$value, value = est$value,
         prior = p0, converged = (est$convergence == 0L),
         convergence = est$convergence, method = "dual",
         call = match.call()),
    class = c("markov_ce", "infometrics")
  )
}

## ---- S3 methods ------------------------------------------------------------
#' @rdname markov_ce
#' @param object,x A \code{markov_ce} object.
#' @param ... Unused.
#' @export
coef.markov_ce <- function(object, ...) object$p
#' @rdname markov_ce
#' @export
fitted.markov_ce <- function(object, ...) object$fitted.values
#' @rdname markov_ce
#' @export
fano_bounds.markov_ce <- function(object, ...)
  .fano_row_bounds(object$p_hat, ncol(object$p_hat))
#' @rdname markov_ce
#' @param average Logical; if \code{TRUE} (default) return the origin-weighted
#'   average marginal effect (\code{S x K}); if \code{FALSE} the full
#'   \code{K x K x S} array of per-origin effects.
#' @param se Logical (default \code{FALSE}); if \code{TRUE} (with
#'   \code{average = TRUE}) attach standard errors.
#' @param se_method \code{"sandwich"} (default, robust and analytic) or
#'   \code{"bootstrap"} (resample units, re-solve). See \code{\link{margins.markov_gce}}.
#' @param B Integer bootstrap resamples when \code{se_method = "bootstrap"}
#'   (default 500).
#' @export
margins.markov_ce <- function(object, average = TRUE, se = FALSE,
                              se_method = c("sandwich", "bootstrap"), B = 500L, ...)
  .margins_markov(object, average, se, match.arg(se_method), B)
#' @rdname markov_ce
#' @export
residuals.markov_ce <- function(object, ...) object$residuals
#' @rdname markov_ce
#' @param digits Significant digits to print.
#' @export
print.markov_ce <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Markov ", x$state_type, " matrix via Cross-Entropy (", length(x$states),
      " states, ", x$n_transitions, " transitions)\n", sep = "")
  cat("\nP (rows sum to 1, P[k, j] = k -> j):\n"); print.default(round(x$p, digits))
  cat("\nmean row entropy H =", format(x$entropy, digits = digits),
      "  normalized S =", format(x$S, digits = digits),
      "  max|moment resid| =", format(x$moment_residual, digits = digits), "\n")
  cat("convergence:", x$convergence,
      if (x$convergence == 0L) "(converged)" else "(see ?optim)", "\n")
  invisible(x)
}
#' @rdname markov_ce
#' @export
summary.markov_ce <- function(object, ...) {
  rows <- apply(object$p, 1L, function(z) { z <- z[z > 0]; -sum(z * log(z)) })
  structure(list(call = object$call, p = object$p, se = object$se, lambda = object$lambda,
                 se_lambda = object$se_lambda, lambda_ref = object$lambda_ref,
                 n_from = object$n_from, states = object$states,
                 n_transitions = object$n_transitions, state_type = object$state_type,
                 row_S = rows / log(length(object$states)),
                 entropy = object$entropy, S = object$S,
                 moment_residual = object$moment_residual,
                 convergence = object$convergence),
            class = "summary.markov_ce")
}
#' @export
print.summary.markov_ce <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Transition matrix P (", x$state_type, ", rows sum to 1), ",
      x$n_transitions, " transitions:\n", sep = "")
  print.default(round(x$p, digits))
  if (x$state_type == "indicator") {
    cat("\nStandard errors (large-sample multinomial):\n"); print.default(round(x$se, digits))
  }
  if (!is.null(rownames(x$lambda))) {
    cat("\nCovariate effects lambda (S x K, original scale):\n"); print.default(round(x$lambda, digits))
  }
  if (!is.null(x$se_lambda) && any(is.finite(x$se_lambda))) {
    cat(sprintf("\nSE(lambda) [analytic-Hessian sandwich, ref destination '%s'; NA = ref / boundary]:\n",
                x$states[x$lambda_ref]))
    print.default(round(x$se_lambda, digits))
  }
  cat("\nOrigin-state weight:", paste(format(x$n_from, digits = digits), collapse = "  "), "\n")
  cat("Per-row normalized entropy S:", paste(format(x$row_S, digits = digits), collapse = "  "), "\n")
  fb <- attr(.fano_row_bounds(x$p, ncol(x$p)), "overall")
  cat("Fano (sec 7.5): mean next-state error pe =", format(fb["mean_pe"], digits = digits),
      ">= bound", format(fb["mean_pe_lower"], digits = digits), "\n")
  cat("Mean row entropy H =", format(x$entropy, digits = digits),
      "  normalized S(P) =", format(x$S, digits = digits),
      "  max|moment resid| =", format(x$moment_residual, digits = digits), "\n")
  cat("convergence:", x$convergence,
      if (x$convergence == 0L) "(converged)" else "(see ?optim)", "\n")
  invisible(x)
}
