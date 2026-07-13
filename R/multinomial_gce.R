# multinomial_gce.R
# nu-weighted Generalized Cross-Entropy multinomial estimator
# (Golan, Judge & Perloff 1996, JASA 91(434):841-853), matrix interface.
#
# Recovers an N x J signal probability matrix p and N x J x M noise weights w
# from a one-hot/share response y (N x J) and design X (N x K), via the
# concentrated dual over the K*(J-1) free Lagrange multipliers (one alternative
# normalized to lambda = 0). Extends the paper with a signal/noise weight nu,
# priors p0/w0, and a which_alternative reference option.
#
#   p_ij  prop to p0_ij  exp(x_i' lambda_j / nu)            (signal)
#   w_ijm prop to w0_ijm exp(x_i' lambda_j v_m / (1-nu))    (noise)
#   e_ij  = sum_m v_m w_ijm
#   dual (MAXIMISED): sum(Y*eta) - nu sum_i logOmega_i - (1-nu) sum_ij logPsi_ij
#   gradient (free block): X'(Y - P - E)[, -which_alternative]
#
# The GJP multinomial GME/GCE estimator (through v0.2.0 it had a formula-
# interface twin, gme_mnl(), since removed); its own argument names per user
# instruction, and the maximise form (fnscale = -1) by user choice. coef()
# returns lambda (the multipliers), not a transformed beta (user instruction).
#
# References:
#   Golan, A., Judge, G. and Perloff, J.M. (1996). A maximum entropy approach
#     to recovering information from multinomial response data. JASA, 91(434),
#     841-853.

# ---- internal: normalization helper ----------------------------------------

#' @keywords internal
.insert_alternative <- function(lambda, which_alternative, kk, jj) {
  z <- matrix(0, nrow = kk, ncol = 1L)              # zero reference column
  if (which_alternative == 1L) {
    cbind(z, lambda)
  } else if (which_alternative < jj) {
    cbind(lambda[, 1:(which_alternative - 1L), drop = FALSE], z,
          lambda[, which_alternative:(jj - 1L), drop = FALSE])
  } else {
    cbind(lambda, z)
  }
}

# ---- internal: default error support --------------------------------------

#' @keywords internal
.mnl_default_support <- function(N, M = 3L) {
  scale <- 1 / sqrt(N)
  seq(-scale, scale, length.out = M)
}

# ---- internal: stable signal / noise pieces --------------------------------
# Returns p (N x J), e (N x J), logOmega (N), logPsi (N x J) for a full K x J
# lambda, using log-sum-exp for numerical stability.

#' @keywords internal
.multinomial_gce_pieces <- function(Lambda, X, p0, w0, v, nu) {
  M1 <- X %*% Lambda                                 # N x J
  ii <- nrow(M1); jj <- ncol(M1); mm <- length(v)

  # signal: p_ij prop to p0_ij exp(M1_ij/nu)  (row-stabilised softmax)
  a    <- log(p0) + M1 / nu                          # N x J
  amax <- apply(a, 1L, max)
  ea   <- exp(a - amax)
  rs   <- rowSums(ea)
  p    <- ea / rs
  logOmega <- amax + log(rs)                         # N

  # noise: w_ijm prop to w0_ijm exp(M1_ij v_m/(1-nu)); per category j
  e      <- matrix(0.0, ii, jj)
  logPsi <- matrix(0.0, ii, jj)
  logw0  <- log(w0)
  for (j in seq_len(jj)) {
    Aj   <- outer(M1[, j], v) / (1 - nu) + logw0[, j, ]   # N x M
    mj   <- apply(Aj, 1L, max)
    eAj  <- exp(Aj - mj)
    rj   <- rowSums(eAj)
    logPsi[, j] <- mj + log(rj)
    e[, j]      <- (eAj / rj) %*% v
  }
  list(p = p, e = e, logOmega = logOmega, logPsi = logPsi)
}

# ---- internal: dual objective and gradient (free params) -------------------

#' @keywords internal
.multinomial_gce_obj <- function(lambda_free, which_alternative, kk, ii, jj,
                                 y, X, p0, w0, v, nu) {
  Lambda <- .insert_alternative(matrix(lambda_free, nrow = kk),
                                which_alternative, kk, jj)
  M1 <- X %*% Lambda
  pc <- .multinomial_gce_pieces(Lambda, X, p0, w0, v, nu)
  sum(y * M1) - nu * sum(pc$logOmega) - (1 - nu) * sum(pc$logPsi)
}

#' @keywords internal
.multinomial_gce_grad <- function(lambda_free, which_alternative, kk, ii, jj,
                                  y, X, p0, w0, v, nu) {
  Lambda <- .insert_alternative(matrix(lambda_free, nrow = kk),
                                which_alternative, kk, jj)
  pc <- .multinomial_gce_pieces(Lambda, X, p0, w0, v, nu)
  g  <- crossprod(X, y - pc$p - pc$e)                # K x J : X'(Y - P - E)
  as.vector(g[, -which_alternative, drop = FALSE])   # drop normalized block
}

# (Fano row-distribution bounds live in the shared .fano_row_bounds() in utils.R.)

# ---- internal: vectorized average marginal effects as a function of the free
# multipliers (for the delta-method Jacobian). Reuses .insert_alternative and the
# row-softmax; returns as.vector(K x J average ME).

#' @keywords internal
.multinomial_gce_me_vec <- function(lambda_free, X, p0, nu, which_alternative,
                                    kk, jj) {
  L    <- .insert_alternative(matrix(lambda_free, nrow = kk),
                              which_alternative, kk, jj)
  a    <- log(p0) + (X %*% L) / nu
  a    <- a - apply(a, 1L, max)
  pp   <- exp(a); pp <- pp / rowSums(pp)
  bar  <- pp %*% t(L)                                      # N x K
  ME   <- matrix(0.0, kk, jj)
  for (j in seq_len(jj))
    ME[, j] <- colMeans((pp[, j] / nu) * sweep(-bar, 2L, L[, j], "+"))
  as.vector(ME)
}

# ---- internal: robust sandwich covariance of the free multipliers -----------
# (-H)^{-1} Vhat (-H)^{-1} = vcov %*% Vhat %*% vcov (since vcov = solve(-H)), with
# Vhat = sum_i g_i g_i' the per-observation dual scores g_i = the free-column block
# of outer(X_i, Y_i - p_i - e_i). Matches the bootstrap; the naive solve(-H)
# delta overstates ~1.35x. Full rank / PD here (N >> K(J-1); no per-obs block).

#' @keywords internal
.multinomial_gce_sandwich <- function(object) {
  X <- object$X; Y <- object$y_mat; p <- object$p_hat; e <- object$e
  ii <- object$N; jj <- object$J; wa <- object$which_alternative
  G  <- matrix(0, ii, nrow(object$vcov))
  for (i in seq_len(ii)) {
    gi <- outer(X[i, ], (Y[i, ] - p[i, ] - e[i, ]))       # K x J
    G[i, ] <- as.vector(gi[, -wa, drop = FALSE])          # free columns
  }
  object$vcov %*% crossprod(G) %*% object$vcov
}

# ===========================================================================
#  multinomial_gce()
# ===========================================================================

#' nu-Weighted Generalized Cross-Entropy Multinomial Estimator (matrix interface)
#'
#' @description
#' Fits an unordered multinomial response model by the Generalized Cross-Entropy
#' estimator of Golan, Judge & Perloff (1996), extended with a signal/noise
#' entropy weight \code{nu}. Recovers an N x J signal probability matrix and the
#' N x J x M noise weights from a one-hot / share response matrix \code{y} and a
#' design matrix \code{X}, via the concentrated dual over the \eqn{K(J-1)} free
#' Lagrange multipliers (one alternative normalized to zero).
#'
#' @details
#' \eqn{p_{ij} \propto p_{0,ij}\exp(x_i'\lambda_j/\nu)} and
#' \eqn{w_{ijm} \propto w_{0,ijm}\exp(x_i'\lambda_j v_m/(1-\nu))}, with
#' \eqn{e_{ij} = \sum_m v_m w_{ijm}}. The dual is \strong{maximized}
#' (\code{fnscale = -1}): \eqn{\sum Y\eta - \nu\sum_i\log\Omega_i -
#' (1-\nu)\sum_{ij}\log\Psi_{ij}}, \eqn{\eta = X\Lambda}, and the free gradient
#' is \eqn{X'(Y - P - E)} dropping the normalized column.
#'
#' At \eqn{\nu = 0.5}, \code{which_alternative = 1}, and a matched support this
#' recovers the same Golan-Judge-Perloff (1996) multinomial GME estimate. The
#' new pieces are \code{nu}, the matrix \code{(y, X)} interface, and
#' \code{which_alternative}. \code{coef()} returns the multipliers
#' \code{lambda} (not a transformed beta).
#'
#' Because the GCE noise forces the reference alternative's error to zero, the
#' arbitrary \code{which_alternative} choice slightly shifts the fitted
#' probabilities (a model property, not a numerical artefact).
#'
#' @param y N x J numeric response matrix: rows are one-hot indicators or
#'   compositional shares summing to 1.
#' @param X N x K design matrix.
#' @param which_alternative Integer in \code{1:J}: the alternative normalized to
#'   \eqn{\lambda = 0} (default 1).
#' @param p0 Optional N x J signal prior (rows sum to 1, strictly positive);
#'   default uniform \code{1/J}.
#' @param w0 Optional N x J x M noise prior (\eqn{(i,j)} fibers sum to 1,
#'   strictly positive); default uniform \code{1/M}.
#' @param v Error support: \code{NULL} (default, a 3-point data-scaled grid on
#'   \eqn{[-1/\sqrt{N}, 1/\sqrt{N}]}), a single whole number giving the support
#'   count, or an explicit numeric vector in \eqn{[-1, 1]}.
#' @param nu Signal/noise entropy weight in \eqn{(0, 1)}; default 0.5.
#' @param control Named list merged over defaults and passed to
#'   \code{\link[stats]{optim}}: \code{maxit} (default 500), \code{reltol}
#'   (default 1e-10). \code{fnscale} is forced to -1 (the dual is maximized).
#'
#' @return An object of class \code{c("multinomial_gce", "infometrics")}. Key
#'   components: \code{lambda}/\code{lambda_hat} (K x J, reference column 0;
#'   \code{coef()} returns it), \code{p}/\code{p_hat}/\code{fitted.values}
#'   (N x J), \code{w}/\code{w_hat} (N x J x M), \code{ee}/\code{e} (N x J),
#'   \code{se}/\code{vcov}, \code{p0}, \code{w0}, \code{v}, \code{nu},
#'   \code{which_alternative}, \code{H_signal} (N-vector), \code{S_p}/\code{S_w}/
#'   \code{S}, \code{value}/\code{objective}, \code{convergence}/\code{converged},
#'   \code{method}, \code{X}, \code{y_mat}, \code{n_misses}, \code{N}/\code{J}/
#'   \code{K}/\code{M}, \code{call}.
#'
#' @references Golan, A., Judge, G. and Perloff, J.M. (1996). A maximum entropy
#'   approach to recovering information from multinomial response data.
#'   \emph{Journal of the American Statistical Association}, \strong{91}(434),
#'   841-853.
#'
#' @seealso \code{\link{margins}} for marginal effects, \code{\link{fano_bounds}}
#'   for information-theoretic error bounds.
#'
#' @examples
#' set.seed(123)
#' n <- 200L; J <- 3L; x <- rnorm(n); X <- cbind(1, x)
#' bt <- cbind(c(0, 0), c(1, 2), c(-0.5, -1))
#' P  <- exp(X %*% bt); P <- P / rowSums(P)
#' yc <- apply(P, 1, function(p) sample(J, 1, prob = p))
#' Y  <- matrix(0, n, J); Y[cbind(seq_len(n), yc)] <- 1
#' fit <- multinomial_gce(Y, X, which_alternative = 1L, nu = 0.5)
#' coef(fit)        # lambda (K x J, reference column 0)
#' margins(fit)              # K x J average marginal effects
#' margins(fit, se = TRUE)   # ... with delta-method standard errors
#' head(fano_bounds(fit))    # Fano error bounds for p_hat (Golan sec 7.5)
#'
#' @importFrom stats optim sd pnorm printCoefmat
#' @export
multinomial_gce <- function(y, X, which_alternative = 1L, p0 = NULL, w0 = NULL,
                            v = NULL, nu = 0.5, control = list()) {
  y <- as.matrix(y); X <- as.matrix(X)
  ii <- nrow(y); jj <- ncol(y); kk <- ncol(X)

  if (nrow(X) != ii)
    stop("nrow(X) (", nrow(X), ") must equal nrow(y) (", ii, ").")
  if (jj < 2L) stop("y must have at least 2 columns (alternatives).")
  if (any(!is.finite(y)) || any(y < 0) || any(abs(rowSums(y) - 1) > 1e-6))
    stop("Rows of 'y' must be non-negative and sum to 1 (one-hot or shares).")
  which_alternative <- as.integer(which_alternative)
  if (which_alternative < 1L || which_alternative > jj)
    stop("which_alternative must be between 1 and ", jj, ".")
  if (!is.numeric(nu) || length(nu) != 1L || nu <= 0 || nu >= 1)
    stop("'nu' must be a single number in (0, 1). Got nu = ", nu, ".")

  # ---- signal prior p0 -----------------------------------------------------
  if (is.null(p0)) p0 <- matrix(1 / jj, ii, jj)
  else {
    p0 <- as.matrix(p0)
    if (!all(dim(p0) == c(ii, jj))) stop("dim(p0) must be ", ii, " x ", jj, ".")
    if (any(p0 <= 0)) stop("All signal prior masses 'p0' must be strictly positive.")
    if (any(abs(rowSums(p0) - 1) > 1e-8)) stop("Rows of 'p0' must sum to 1.")
  }

  # ---- error support v (scalar = count; vector = support in [-1,1]) --------
  if (is.null(v)) {
    mm <- if (is.null(w0)) 3L else dim(w0)[3]
    v  <- .mnl_default_support(ii, mm)              # 3-point grid on +/- 1/sqrt(N)
  } else if (length(v) == 1L) {
    mm <- as.integer(v); if (mm < 2L) stop("scalar 'v' (support count) must be >= 2.")
    v  <- .mnl_default_support(ii, mm)
  } else {
    v <- as.numeric(v)
    if (any(v < -1 | v > 1)) stop("All elements of 'v' must lie within [-1, 1].")
    mm <- length(v)
  }

  # ---- noise prior w0 ------------------------------------------------------
  if (is.null(w0)) w0 <- array(1 / mm, dim = c(ii, jj, mm))
  else {
    w0 <- as.array(w0)
    if (!all(dim(w0) == c(ii, jj, mm)))
      stop("dim(w0) must be ", ii, " x ", jj, " x ", mm, ".")
    if (any(w0 <= 0)) stop("All noise prior masses 'w0' must be strictly positive.")
    if (any(abs(apply(w0, c(1, 2), sum) - 1) > 1e-8))
      stop("Each (i,j) fiber of w0 must sum to 1.")
  }

  con <- list(maxit = 500L, reltol = 1e-10)
  con[names(control)] <- control
  con$fnscale <- -1                                 # maximise the dual

  # ---- optimise ------------------------------------------------------------
  lambda0 <- rep(0.0, kk * (jj - 1L))
  est <- stats::optim(
    par = lambda0, fn = .multinomial_gce_obj, gr = .multinomial_gce_grad,
    which_alternative = which_alternative, kk = kk, ii = ii, jj = jj,
    y = y, X = X, p0 = p0, w0 = w0, v = v, nu = nu,
    method = "BFGS", control = con, hessian = TRUE
  )
  if (est$convergence != 0L)
    warning("optim did not converge (code ", est$convergence, ").")

  Lambda <- .insert_alternative(matrix(est$par, nrow = kk),
                                which_alternative, kk, jj)
  pc <- .multinomial_gce_pieces(Lambda, X, p0, w0, v, nu)
  p  <- pc$p; ee <- pc$e

  # ---- SE of lambda (maximise => -Hessian is the information matrix) --------
  vcov_free <- tryCatch(
    solve(-est$hessian),
    error = function(e) {
      warning("Hessian singular; SEs set to NA.")
      matrix(NA_real_, nrow(est$hessian), ncol(est$hessian))
    }
  )
  se_free   <- sqrt(pmax(diag(vcov_free), 0))
  SE        <- matrix(0.0, kk, jj)
  SE[, -which_alternative] <- matrix(se_free, kk, jj - 1L)

  # ---- noise weights (N x J x M) for w_hat ---------------------------------
  W <- array(0.0, dim = c(ii, jj, mm)); logw0 <- log(w0)
  M1 <- X %*% Lambda
  for (j in seq_len(jj)) {
    Aj <- outer(M1[, j], v) / (1 - nu) + logw0[, j, ]
    mj <- apply(Aj, 1L, max); eAj <- exp(Aj - mj)
    W[, j, ] <- eAj / rowSums(eAj)
  }

  # ---- names ---------------------------------------------------------------
  rn <- colnames(X); cn <- colnames(y)
  if (is.null(rn)) rn <- paste0("x", seq_len(kk))
  if (is.null(cn)) cn <- paste0("alt", seq_len(jj))
  dimnames(Lambda) <- dimnames(SE) <- list(rn, cn)
  dimnames(p) <- dimnames(ee) <- list(NULL, cn)

  # ---- entropy diagnostics -------------------------------------------------
  Plog <- ifelse(p > 0, p * log(p), 0.0)
  H_signal <- -rowSums(Plog)                        # N-vector
  S_p  <- -sum(Plog) / (log(jj) * ii)
  Wf   <- as.numeric(W); Wlog <- ifelse(Wf > 0, Wf * log(Wf), 0.0)
  S_w  <- -sum(Wlog) / (log(mm) * ii * jj)
  n_misses <- sum(max.col(p) != max.col(y))

  structure(
    list(
      # draft names
      lambda = Lambda, p = p, w = W, ee = ee,
      p0 = p0, w0 = w0, v = v, nu = nu,
      value = est$value, convergence = est$convergence,
      # CLAUDE.md canonical aliases
      lambda_hat = Lambda, p_hat = p, w_hat = W, e = ee,
      objective = est$value, converged = (est$convergence == 0L),
      method = "dual", H_signal = H_signal, S = S_p + S_w,
      # additions
      fitted.values = p, y_mat = y, X = X,
      se = SE, vcov = vcov_free, which_alternative = which_alternative,
      S_p = S_p, S_w = S_w, n_misses = n_misses, hessian = est$hessian,
      N = ii, J = jj, K = kk, M = mm, call = match.call()
    ),
    class = c("multinomial_gce", "infometrics")
  )
}

# ===========================================================================
#  margins() generic + method
# ===========================================================================

#' Marginal Effects
#'
#' Generic for average (or per-observation) marginal effects of a fitted model.
#'
#' @param object A fitted model object.
#' @param ... Passed to methods.
#' @return A method-specific object of marginal effects.
#' @seealso \code{\link{margins.multinomial_gce}}
#' @export
margins <- function(object, ...) UseMethod("margins")

#' Marginal Effects for a Multinomial GCE Fit
#'
#' Partial derivatives of the signal probabilities with respect to each column
#' of \code{X}:
#' \deqn{\partial p_{ij}/\partial x_{ik}
#'   = (p_{ij}/\nu)\,(\lambda_{kj} - \textstyle\sum_l p_{il}\lambda_{kl}).}
#' The reference alternative's column is included like any other; its
#' "intercept" effect (derivative w.r.t. a constant) is rarely interpreted.
#'
#' @details
#' With \code{se = TRUE} (and \code{average = TRUE}) standard errors for the
#' average marginal effects are attached. These are a \emph{sampling} quantity,
#' \strong{unrelated} to the Fano \code{\link{fano_bounds}} (which bound
#' classification error, not sampling variance). All methods sandwich a numerical
#' Jacobian \eqn{J} of the average effects in the free multipliers with a
#' covariance of \eqn{\hat\lambda}:
#' \itemize{
#'   \item \code{se_method = "sandwich"} (default) uses the robust covariance
#'     \eqn{(-H)^{-1}\hat V(-H)^{-1}} with \eqn{\hat V=\sum_i g_i g_i'} the
#'     per-observation dual scores. This matches a row bootstrap and is the
#'     accurate choice under the GME's finite-support regularization.
#'   \item \code{se_method = "delta"} uses the naive Hessian inverse
#'     \eqn{V=}\code{vcov(object)} (the classical Golan 2008 sec. 7.3
#'     observed-information SE); being tied to \eqn{(-H)^{-1}} it is
#'     \emph{conservative} (~1.35x the bootstrap here).
#'   \item \code{se_method = "bootstrap"} resamples observations, refits, and
#'     takes the across-resample SD (accurate, slower).
#' }
#'
#' @param object A \code{multinomial_gce} object.
#' @param average Logical (default \code{TRUE}). If \code{TRUE}, return the
#'   K x J matrix of marginal effects averaged over observations; if
#'   \code{FALSE}, the full N x J x K array of per-observation derivatives.
#' @param se Logical (default \code{FALSE}). If \code{TRUE} (requires
#'   \code{average = TRUE}), return a \code{margins_gce} object carrying the
#'   estimates, standard errors, z-values, and p-values.
#' @param se_method \code{"sandwich"} (default), \code{"delta"}, or
#'   \code{"bootstrap"}; see Details.
#' @param B Integer number of bootstrap resamples when
#'   \code{se_method = "bootstrap"} (default 500).
#' @param ... Unused.
#' @return With \code{se = FALSE}: a K x J matrix (average effects) or an
#'   N x J x K array. With \code{se = TRUE}: a \code{margins_gce} object.
#' @seealso \code{\link{fano_bounds}} for information-theoretic error bounds.
#' @export
margins.multinomial_gce <- function(object, average = TRUE, se = FALSE,
                                    se_method = c("sandwich", "delta", "bootstrap"),
                                    B = 500L, ...) {
  p <- object$p_hat; lambda <- object$lambda; nu <- object$nu
  N <- nrow(p); J <- ncol(p); K <- nrow(lambda)
  bar <- p %*% t(lambda)                            # N x K: sum_l p_il lambda[k,l]

  me_average <- function() {
    out <- matrix(0.0, K, J,
                  dimnames = list(rownames(lambda), colnames(lambda)))
    for (j in seq_len(J))
      out[, j] <- colMeans((p[, j] / nu) * sweep(-bar, 2L, lambda[, j], "+"))
    out
  }

  if (!se) {
    if (average) return(me_average())
    ME <- array(0.0, dim = c(N, J, K),
                dimnames = list(NULL, colnames(lambda), rownames(lambda)))
    for (j in seq_len(J))
      ME[, j, ] <- (p[, j] / nu) * sweep(-bar, 2L, lambda[, j], "+")
    return(ME)
  }

  # ---- se = TRUE: average marginal effects with standard errors -------------
  if (!average)
    stop("Standard errors are implemented only for the average effects ",
         "(average = TRUE).")
  se_method <- match.arg(se_method)
  wa   <- object$which_alternative
  ME   <- me_average()

  if (se_method != "bootstrap") {
    # sandwich (default) and delta share the ME Jacobian; only the covariance of
    # the free multipliers differs.
    Vlam <- if (se_method == "sandwich") .multinomial_gce_sandwich(object)
            else object$vcov                        # "delta": naive Hessian inverse
    lf  <- as.vector(lambda[, -wa, drop = FALSE])   # free multipliers
    h   <- 1e-5; nf <- length(lf)
    Jac <- matrix(0.0, K * J, nf)
    for (m in seq_len(nf)) {
      e <- numeric(nf); e[m] <- h
      Jac[, m] <- (.multinomial_gce_me_vec(lf + e, object$X, object$p0, nu, wa, K, J) -
                   .multinomial_gce_me_vec(lf - e, object$X, object$p0, nu, wa, K, J)) / (2 * h)
    }
    Vme <- Jac %*% Vlam %*% t(Jac)
    SE  <- matrix(sqrt(pmax(diag(Vme), 0)), K, J, dimnames = dimnames(lambda))
  } else {
    MEb <- array(NA_real_, dim = c(K, J, B))
    for (b in seq_len(B)) {
      idx <- sample(N, N, replace = TRUE)
      fb  <- tryCatch(
        multinomial_gce(object$y_mat[idx, , drop = FALSE],
                        object$X[idx, , drop = FALSE],
                        which_alternative = wa,
                        p0 = object$p0[idx, , drop = FALSE],
                        w0 = object$w0[idx, , , drop = FALSE],
                        v  = object$v, nu = nu),
        error = function(e) NULL)
      if (is.null(fb) || !fb$converged) next
      MEb[, , b] <- matrix(
        .multinomial_gce_me_vec(as.vector(fb$lambda[, -wa, drop = FALSE]),
                                object$X[idx, , drop = FALSE], fb$p0, nu, wa, K, J),
        K, J)
    }
    n_ok <- sum(!is.na(MEb[1, 1, ]))
    if (n_ok < 2L) stop("bootstrap produced < 2 usable resamples.")
    SE <- matrix(apply(MEb, c(1L, 2L), stats::sd, na.rm = TRUE), K, J,
                 dimnames = dimnames(lambda))
    attr(SE, "n_boot") <- n_ok
  }

  z <- ME / SE
  structure(
    list(estimate = ME, se = SE, z = z, p.value = 2 * pnorm(-abs(z)),
         se_method = se_method, which_alternative = wa,
         n_boot = attr(SE, "n_boot")),
    class = "margins_gce"
  )
}

#' @export
print.margins_gce <- function(x, digits = 4L, ...) {
  cat(sprintf("\nAverage marginal effects (SE: %s)\n", x$se_method))
  if (!is.null(x$n_boot))
    cat(sprintf("  bootstrap resamples used: %d\n", x$n_boot))
  cat(strrep("-", 48L), "\n")
  cn <- colnames(x$estimate)
  for (j in seq_len(ncol(x$estimate))) {
    cat(sprintf("\n  d p(%s) / d x:\n", cn[j]))
    tab <- cbind(Estimate = x$estimate[, j], `Std. Error` = x$se[, j],
                 `z value` = x$z[, j], `Pr(>|z|)` = x$p.value[, j])
    printCoefmat(tab, digits = digits, has.Pvalue = TRUE)
  }
  invisible(x)
}

# ===========================================================================
#  fano_bounds() generic + method  (Golan 2008, sec 7.5)
# ===========================================================================

#' Information-Theoretic (Fano) Error Bounds
#'
#' Generic for Fano's-inequality error bounds on a fitted probability matrix.
#'
#' @param object A fitted model object.
#' @param ... Passed to methods.
#' @return A method-specific object of error bounds.
#' @seealso \code{\link{fano_bounds.multinomial_gce}}
#' @export
fano_bounds <- function(object, ...) UseMethod("fano_bounds")

#' Fano Error Bounds for a Multinomial GCE Fit
#'
#' For each observation (a row of \code{p_hat}, a distribution over the \eqn{J}
#' alternatives) the modal classifier predicts \eqn{\arg\max_j p_{ij}} with error
#' probability \eqn{pe_i = 1 - \max_j p_{ij}}. Fano's inequality (Golan 2008,
#' sec 7.5) lower-bounds this error by the normalized entropy:
#' \deqn{pe_i \ge S(p_i) - \log 2/\log J,\qquad S(p_i)=H(p_i)/\log J.}
#' This is an \strong{information-theoretic error bound} on prediction accuracy,
#' \emph{not} a sampling standard error (for the latter see
#' \code{\link{margins}} with \code{se = TRUE}).
#'
#' @param object A \code{multinomial_gce} object.
#' @param ... Unused.
#' @return A data frame with one row per observation and columns \code{p_max},
#'   \code{pe} (modal error), \code{H} (entropy, nats), \code{S} (normalized
#'   entropy), and \code{pe_lower} (Fano weak lower bound). An \code{"overall"}
#'   attribute holds \code{mean_pe}, \code{mean_pe_lower}, and \code{S_system}.
#' @references Golan, A. (2008). Information and Entropy Econometrics.
#'   \emph{Foundations and Trends in Econometrics}, \strong{2}(1-2), 1-145
#'   (sec 3.6, 7.5); Fano, R. (1961). \emph{Transmission of Information}.
#' @seealso \code{\link{margins}}
#' @export
fano_bounds.multinomial_gce <- function(object, ...) {
  .fano_row_bounds(object$p_hat, object$J)
}

# ===========================================================================
#  S3 methods
# ===========================================================================

#' @export
coef.multinomial_gce <- function(object, ...) object$lambda

#' @export
fitted.multinomial_gce <- function(object, ...) object$p_hat

#' @export
residuals.multinomial_gce <- function(object, ...)
  object$y_mat - object$p_hat - object$e

#' @export
print.multinomial_gce <- function(x, digits = 4L, ...) {
  uniform_p <- max(abs(x$p0 - 1 / x$J)) < 1e-10
  est <- if (uniform_p) "GME" else "GCE"
  cat(sprintf("\nMultinomial %s (nu-weighted, GJP 1996)\n", est))
  cat(strrep("-", 48L), "\n")
  cat(sprintf("  N = %d   J = %d   K = %d   M = %d   nu = %.3g\n",
              x$N, x$J, x$K, x$M, x$nu))
  cat(sprintf("  support v  : [%.4f, %.4f]   reference = alt %d\n",
              min(x$v), max(x$v), x$which_alternative))
  cat(sprintf("  S(p)=%.4f  S(w)=%.4f  misses=%d/%d  convergence=%s\n",
              x$S_p, x$S_w, x$n_misses, x$N,
              if (x$converged) "yes" else "no"))
  cat("\nlambda (reference column = 0):\n")
  print(round(x$lambda, digits))
  invisible(x)
}

#' @export
summary.multinomial_gce <- function(object, digits = 4L, ...) {
  J <- object$J; ref <- object$which_alternative
  cols <- setdiff(seq_len(J), ref)
  tabs <- lapply(cols, function(j) data.frame(
    lambda = object$lambda[, j],
    se     = object$se[, j],
    t      = ifelse(object$se[, j] > 0, object$lambda[, j] / object$se[, j], NA_real_)
  ))
  names(tabs) <- colnames(object$lambda)[cols]

  cat(sprintf("\nMultinomial GCE (nu = %.3g, GJP 1996)\n", object$nu))
  cat("Call: "); print(object$call); cat("\n")
  cat(sprintf("  N=%d  J=%d  K=%d  M=%d  reference=alt %d\n",
              object$N, object$J, object$K, object$M, ref))
  cat(sprintf("  S(p)=%.4f  S(w)=%.4f  misses=%d/%d  convergence=%s\n",
              object$S_p, object$S_w, object$n_misses, object$N,
              if (object$converged) "yes" else "no"))
  fb <- attr(.fano_row_bounds(object$p_hat, J), "overall")
  cat(sprintf("  Fano (sec 7.5): mean modal error pe=%.4f  >= mean bound %.4f  [S(P)=%.4f]\n",
              fb["mean_pe"], fb["mean_pe_lower"], fb["S_system"]))
  for (nm in names(tabs)) {
    cat(sprintf("\n  Alternative '%s' (vs. reference '%s'):\n",
                nm, colnames(object$lambda)[ref]))
    print(round(tabs[[nm]], digits))
  }
  invisible(object)
}
