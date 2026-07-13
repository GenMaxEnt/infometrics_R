# mixed_gce.R
# Doubly-reparameterized Generalized Cross-Entropy "mixed" model.
#
# Response Y (N x J shares), design X (N x J x K). The signal probability is
# itself reparameterized on a support s in [0,1]:
#   theta_ijm prop to theta0_ijm exp( s_m (V_ij + rho_i) / nu ),  V_ij = sum_k X_ijk lambda_kj
#   p_ij = sum_m s_m theta_ijm                       (signal, with sum_j p_ij = 1 via rho_i)
#   w_ijh prop to w0_ijh exp( u_h V_ij / (1-nu) ),   e_ij = sum_h u_h w_ijh   (noise)
# solved via the concentrated dual (MAXIMISED, fnscale = -1) over (rho, lambda):
#   LL = sum_ij Y_ij V_ij + sum_i rho_i
#        - nu sum_ij log Omega_ij - (1-nu) sum_ij log Psi_ij
#   grad_rho_i  = 1 - sum_j p_ij                      (adding-up constraint)
#   grad_lam_kj = sum_i X_ijk (Y_ij - p_ij - e_ij)    (data moment)
#
# Genuinely new estimator: distinct from multinomial_gce() (which recovers p
# directly with a reference-column normalization); here p is a support-weighted
# mean of theta and adding-up is a separate rho block. Keeps the draft's own
# argument names (Y, X, theta0, w0, nu) plus exposed supports s, u. coef()
# returns lambda (GJP-family convention), and the maximise form is kept (with
# markov_*, multinomial_gce, linreg_iv, panel_gce).
#
# References:
#   Golan, A., Judge, G. and Perloff, J.M. (1996). A maximum entropy approach to
#     recovering information from multinomial response data. JASA, 91, 841-853.
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145.

# ---- internal: V_ij = sum_k X_ijk lambda_kj --------------------------------

#' @keywords internal
.mixed_gce_V <- function(X, lambda, ii, jj) {
  V <- matrix(0, ii, jj)
  for (j in seq_len(jj)) V[, j] <- matrix(X[, j, ], nrow = ii) %*% lambda[, j]
  V
}

# ---- internal: signal pieces (log-sum-exp stabilised) ----------------------
# theta (N x J x M), p (N x J), logOmega (N x J).

#' @keywords internal
.mixed_gce_signal <- function(XL, s, theta0, nu) {
  d <- dim(theta0); mm <- d[3]
  A <- array(0, d)
  logt0 <- log(theta0)
  for (m in seq_len(mm)) A[, , m] <- logt0[, , m] + s[m] * XL / nu
  amax <- apply(A, c(1L, 2L), max)
  eA   <- exp(sweep(A, c(1L, 2L), amax, "-"))
  rs   <- rowSums(eA, dims = 2L)
  theta <- sweep(eA, c(1L, 2L), rs, "/")
  p    <- rowSums(sweep(theta, 3L, s, "*"), dims = 2L)
  list(theta = theta, p = p, logOmega = amax + log(rs))
}

# ---- internal: noise pieces (log-sum-exp stabilised) -----------------------
# w (N x J x H), e (N x J), logPsi (N x J).

#' @keywords internal
.mixed_gce_noise <- function(V, u, w0, nu) {
  d <- dim(w0); hh <- d[3]
  B <- array(0, d)
  logw0 <- log(w0)
  for (h in seq_len(hh)) B[, , h] <- logw0[, , h] + u[h] * V / (1 - nu)
  bmax <- apply(B, c(1L, 2L), max)
  eB   <- exp(sweep(B, c(1L, 2L), bmax, "-"))
  rs   <- rowSums(eB, dims = 2L)
  w    <- sweep(eB, c(1L, 2L), rs, "/")
  e    <- rowSums(sweep(w, 3L, u, "*"), dims = 2L)
  list(w = w, e = e, logPsi = bmax + log(rs))
}

# ---- internal: dual objective and gradient ---------------------------------

#' @keywords internal
.mixed_gce_obj <- function(param, Y, X, ii, jj, kk, s, u, theta0, w0, nu) {
  rho    <- param[1:ii]
  lambda <- matrix(param[(ii + 1L):(ii + kk * jj)], kk, jj)
  V  <- .mixed_gce_V(X, lambda, ii, jj)
  XL <- V + rho                                      # recycles rho down the rows
  sg <- .mixed_gce_signal(XL, s, theta0, nu)
  ns <- .mixed_gce_noise(V, u, w0, nu)
  sum(Y * V) + sum(rho) - nu * sum(sg$logOmega) - (1 - nu) * sum(ns$logPsi)
}

#' @keywords internal
.mixed_gce_grad <- function(param, Y, X, ii, jj, kk, s, u, theta0, w0, nu) {
  rho    <- param[1:ii]
  lambda <- matrix(param[(ii + 1L):(ii + kk * jj)], kk, jj)
  V  <- .mixed_gce_V(X, lambda, ii, jj)
  XL <- V + rho
  sg <- .mixed_gce_signal(XL, s, theta0, nu)
  ns <- .mixed_gce_noise(V, u, w0, nu)
  g_rho <- 1 - rowSums(sg$p)                         # length ii
  R     <- Y - sg$p - ns$e                           # ii x jj
  g_lam <- matrix(0, kk, jj)
  for (j in seq_len(jj))
    g_lam[, j] <- colSums(matrix(X[, j, ], nrow = ii) * R[, j])   # sum_i X_ijk R_ij
  c(g_rho, as.vector(g_lam))
}

# ===========================================================================
#  mixed_gce()
# ===========================================================================

#' Doubly-Reparameterized Generalized Cross-Entropy Mixed Model
#'
#' @description
#' Fits the Generalized Cross-Entropy "mixed" model in which the signal
#' probabilities are themselves reparameterized on a bounded support. From a
#' share response \code{Y} (N x J) and a design array \code{X} (N x J x K) it
#' recovers the coefficients \eqn{\lambda} (K x J), the per-observation adding-up
#' multipliers \eqn{\rho} (N), the signal weights \eqn{\theta} (N x J x M) with
#' \eqn{p_{ij}=\sum_m s_m\theta_{ijm}}, and the noise weights \eqn{w} (N x J x H)
#' with \eqn{e_{ij}=\sum_h u_h w_{ijh}}, via the concentrated dual.
#'
#' @details
#' With \eqn{V_{ij}=\sum_k X_{ijk}\lambda_{kj}},
#' \eqn{\theta_{ijm}\propto\theta^0_{ijm}\exp[s_m(V_{ij}+\rho_i)/\nu]} and
#' \eqn{w_{ijh}\propto w^0_{ijh}\exp[u_h V_{ij}/(1-\nu)]}. The dual is
#' \strong{maximized} (\code{fnscale = -1}):
#' \deqn{LL = \sum_{ij} Y_{ij}V_{ij} + \sum_i \rho_i
#'   - \nu\sum_{ij}\log\Omega_{ij} - (1-\nu)\sum_{ij}\log\Psi_{ij},}
#' whose first-order conditions are the adding-up constraint
#' \eqn{1-\sum_j p_{ij}=0} (in \eqn{\rho_i}) and the data moment
#' \eqn{\sum_i X_{ijk}(Y_{ij}-p_{ij}-e_{ij})=0} (in \eqn{\lambda_{kj}}). Uniform
#' priors give the GME (maximum-entropy) solution; user priors give GCE. As the
#' noise support \eqn{u\to 0} the fit approaches the pure-signal solution.
#'
#' @param Y N x J numeric response matrix of shares (rows are compositional,
#'   summing to about 1 up to the noise term).
#' @param X N x J x K design array (\code{X[i, j, ]} is the covariate vector for
#'   observation \eqn{i} and category \eqn{j}).
#' @param theta0 Optional N x J x M signal prior (each \eqn{(i,j)} fiber > 0 and
#'   summing to 1); default uniform \code{1/M}. \code{M} is taken from \code{s}.
#' @param w0 Optional N x J x H noise prior (each \eqn{(i,j)} fiber > 0 and
#'   summing to 1); default uniform \code{1/H}. \code{H} is taken from \code{u}.
#' @param s Signal support in \eqn{[0,1]}: \code{NULL} (default 3-point grid
#'   \code{c(0, 0.5, 1)}), a single whole number giving the number of
#'   equally-spaced points, or an explicit numeric vector in \eqn{[0,1]}.
#' @param u Noise support in \eqn{[-1,1]}: \code{NULL} (default 3-point grid
#'   \code{c(-1, 0, 1)}), a single whole number giving the number of
#'   equally-spaced points, or an explicit symmetric numeric vector in
#'   \eqn{[-1,1]}. Supply your own \code{u} to widen or tighten the noise term
#'   independently of \code{s} and the priors.
#' @param nu Signal/noise entropy weight in \eqn{(0,1)}; default 0.5.
#' @param control Named list merged over defaults and passed to
#'   \code{\link[stats]{optim}} (BFGS): \code{maxit} (default 1000) and
#'   \code{reltol} (default 1e-12). \code{fnscale} is forced to -1.
#'
#' @return An object of class \code{c("mixed_gce", "infometrics")} with
#'   \code{lambda}/\code{lambda_hat} (K x J; \code{coef()} returns it),
#'   \code{rho}/\code{rho_hat} (N), \code{theta}/\code{theta_hat} (N x J x M),
#'   \code{p}/\code{p_hat}/\code{fitted.values} (N x J), \code{w}/\code{w_hat}
#'   (N x J x H), \code{e}/\code{e_hat} (N x J), \code{V}, standard errors
#'   \code{se_lambda} (K x J)/\code{se_rho} (N)/\code{vcov} (via the dual
#'   Hessian), the stored inputs and supports, entropy diagnostics
#'   \code{H_signal} (N), \code{S_p}/\code{S_w}/\code{S},
#'   \code{objective}/\code{value}, \code{converged}/\code{convergence},
#'   \code{method}, dims \code{N}/\code{J}/\code{K}/\code{M}/\code{H}, and
#'   \code{call}.
#'
#' @references
#' Golan, A., Judge, G. and Perloff, J.M. (1996). A maximum entropy approach to
#' recovering information from multinomial response data. \emph{Journal of the
#' American Statistical Association}, \strong{91}(434), 841-853.
#'
#' @seealso \code{\link{multinomial_gce}} (the reference-normalized sibling),
#'   \code{\link{margins}} for marginal effects.
#'
#' @examples
#' set.seed(1)
#' N <- 25L; J <- 3L; K <- 2L
#' X <- array(0, dim = c(N, J, K))
#' X[, , 1] <- 1                                   # intercept-like column
#' X[, , 2] <- matrix(rnorm(N * J), N, J)          # a covariate
#' Y <- matrix(runif(N * J), N, J); Y <- Y / rowSums(Y)   # compositional shares
#' fit <- mixed_gce(Y, X, nu = 0.5)
#' coef(fit)                                        # lambda (K x J)
#' margins(fit)                                     # K x J average marginal effects
#' margins(fit, se = TRUE)                          # ... with robust sandwich SEs
#' head(fano_bounds(fit))                           # Fano error bounds (Golan sec 7.5)
#'
#' @importFrom stats optim sd pnorm printCoefmat
#' @export
mixed_gce <- function(Y, X, theta0 = NULL, w0 = NULL, s = NULL, u = NULL,
                      nu = 0.5, control = list()) {
  Y <- as.matrix(Y); X <- as.array(X)
  if (length(dim(X)) != 3L)
    stop("X must be a 3-dimensional array (N x J x K). Got ", length(dim(X)),
         " dimension(s).")
  d <- dim(X); ii <- d[1]; jj <- d[2]; kk <- d[3]
  if (nrow(Y) != ii || ncol(Y) != jj)
    stop("dim(Y) must be ", ii, " x ", jj, " (the first two dims of X). Got ",
         nrow(Y), " x ", ncol(Y), ".")
  if (anyNA(Y) || anyNA(X)) stop("'Y' and 'X' must not contain NA.")
  if (!is.numeric(nu) || length(nu) != 1L || nu <= 0 || nu >= 1)
    stop("'nu' must be a single number strictly inside (0, 1). Got nu = ", nu, ".")

  # ---- signal support s on [0, 1] ------------------------------------------
  if (is.null(s)) {
    s <- c(0, 0.5, 1)
  } else if (length(s) == 1L) {
    mm <- as.integer(s)
    if (mm < 2L) stop("scalar 's' (support count) must be >= 2. Got ", mm, ".")
    s <- seq(0, 1, length.out = mm)
  } else {
    s <- as.numeric(s)
    if (any(s < 0 | s > 1)) stop("All elements of 's' must lie within [0, 1].")
  }
  mm <- length(s)

  # ---- noise support u symmetric on [-1, 1] --------------------------------
  if (is.null(u)) {
    u <- c(-1, 0, 1)
  } else if (length(u) == 1L) {
    hh <- as.integer(u)
    if (hh < 2L) stop("scalar 'u' (support count) must be >= 2. Got ", hh, ".")
    u <- seq(-1, 1, length.out = hh)
  } else {
    u <- as.numeric(u)
    if (any(u < -1 | u > 1)) stop("All elements of 'u' must lie within [-1, 1].")
    if (abs(mean(range(u))) > 1e-8 * max(diff(range(u)), 1))
      stop("Noise support 'u' must be symmetric about zero.")
  }
  hh <- length(u)

  # ---- priors --------------------------------------------------------------
  fix_prior <- function(pr, name, mmhh) {
    if (is.null(pr)) return(array(1 / mmhh, dim = c(ii, jj, mmhh)))
    pr <- as.array(pr)
    if (!all(dim(pr) == c(ii, jj, mmhh)))
      stop("dim(", name, ") must be ", ii, " x ", jj, " x ", mmhh, ".")
    if (any(pr <= 0)) stop("All entries of '", name, "' must be strictly positive.")
    if (any(abs(apply(pr, c(1L, 2L), sum) - 1) > 1e-8))
      stop("Each (i, j) fiber of '", name, "' must sum to 1.")
    pr
  }
  theta0 <- fix_prior(theta0, "theta0", mm)
  w0     <- fix_prior(w0, "w0", hh)

  # ---- optimise the concentrated dual --------------------------------------
  con <- list(maxit = 1000L, reltol = 1e-12)
  con[names(control)] <- control
  con$fnscale <- -1
  param0 <- rep(0, ii + kk * jj)
  est <- stats::optim(param0, .mixed_gce_obj, .mixed_gce_grad,
                      Y = Y, X = X, ii = ii, jj = jj, kk = kk, s = s, u = u,
                      theta0 = theta0, w0 = w0, nu = nu,
                      method = "BFGS", control = con, hessian = TRUE)
  if (est$convergence != 0L)
    warning("optim did not converge (code ", est$convergence, ").")

  rho    <- est$par[1:ii]
  lambda <- matrix(est$par[(ii + 1L):(ii + kk * jj)], kk, jj)
  V  <- .mixed_gce_V(X, lambda, ii, jj); XL <- V + rho
  sg <- .mixed_gce_signal(XL, s, theta0, nu)
  ns <- .mixed_gce_noise(V, u, w0, nu)
  p  <- sg$p; theta <- sg$theta; w <- ns$w; e <- ns$e

  # ---- standard errors: -Hessian of the (maximised) dual is the info matrix -
  vcov <- tryCatch(solve(-est$hessian), error = function(err) {
    warning("dual Hessian singular; standard errors set to NA.")
    matrix(NA_real_, length(est$par), length(est$par))
  })
  se_all    <- sqrt(pmax(diag(vcov), 0))
  se_rho    <- se_all[1:ii]
  se_lambda <- matrix(se_all[(ii + 1L):(ii + kk * jj)], kk, jj)

  # ---- names ---------------------------------------------------------------
  rn <- dimnames(X)[[3]]; cn <- colnames(Y)
  if (is.null(rn)) rn <- paste0("x", seq_len(kk))
  if (is.null(cn)) cn <- paste0("alt", seq_len(jj))
  dimnames(lambda) <- dimnames(se_lambda) <- list(rn, cn)
  dimnames(p) <- list(NULL, cn)

  # ---- entropy diagnostics -------------------------------------------------
  tlog <- ifelse(theta > 0, theta * log(theta), 0)
  H_signal <- apply(tlog, 1L, function(z) -sum(z))   # N-vector (over j, m)
  S_p  <- -sum(tlog) / (ii * jj * log(mm))
  wlog <- ifelse(w > 0, w * log(w), 0)
  S_w  <- -sum(wlog) / (ii * jj * log(hh))

  structure(
    list(
      # draft names
      rho_hat = rho, lambda_hat = lambda, theta_hat = theta,
      p_hat = p, w_hat = w, e_hat = e, convergence = est$convergence,
      value = est$value,
      # canonical aliases
      lambda = lambda, rho = rho, theta = theta, p = p, w = w, e = e,
      fitted.values = p, objective = est$value,
      converged = (est$convergence == 0L), method = "dual",
      H_signal = H_signal, S = S_p, S_p = S_p, S_w = S_w,
      # standard errors
      vcov = vcov, se_lambda = se_lambda, se_rho = se_rho, hessian = est$hessian,
      # stored inputs / dims
      V = V, X = X, y_mat = Y, s = s, u = u, theta0 = theta0, w0 = w0, nu = nu,
      N = ii, J = jj, K = kk, M = mm, H = hh, call = match.call()
    ),
    class = c("mixed_gce", "infometrics")
  )
}

# ===========================================================================
#  margins() method
# ===========================================================================

# ---- internal: vectorized average marginal effects as a function of the full
# parameter vector (rho, lambda), for the sandwich/bootstrap Jacobian. ---------

#' @keywords internal
.mixed_gce_me_vec <- function(param, X, s, theta0, nu, ii, jj, kk) {
  rho    <- param[1:ii]
  lambda <- matrix(param[(ii + 1L):(ii + kk * jj)], kk, jj)
  V   <- .mixed_gce_V(X, lambda, ii, jj)
  sg  <- .mixed_gce_signal(V + rho, s, theta0, nu)
  Es2 <- apply(sweep(sg$theta, 3L, s^2, "*"), c(1L, 2L), sum)
  vbar <- colMeans(Es2 - sg$p^2)                            # J
  ME   <- matrix(0, kk, jj)
  for (j in seq_len(jj)) ME[, j] <- (lambda[, j] / nu) * vbar[j]
  as.vector(ME)
}

# ---- internal: robust sandwich covariance of (rho, lambda) -------------------
# (-H)^{-1} Vhat (-H)^{-1} with Vhat = sum_i g_i g_i' the per-observation dual
# scores: rho-block = (1 - sum_j p_ij) in slot i (~0 at the optimum, so rho is
# pinned by the adding-up constraint); lambda-block = X_ijk (Y_ij - p_ij - e_ij).
# This matches the bootstrap; the naive solve(-H) overstates ~7-9x here.

#' @keywords internal
.mixed_gce_sandwich <- function(object) {
  X <- object$X; Y <- object$y_mat; p <- object$p_hat; e <- object$e_hat
  ii <- object$N; jj <- object$J; kk <- object$K
  np <- ii + kk * jj
  G  <- matrix(0, ii, np)
  for (i in seq_len(ii)) {
    G[i, i] <- 1 - sum(p[i, ])                              # rho_i score
    glam <- matrix(0, kk, jj)
    for (j in seq_len(jj)) glam[, j] <- X[i, j, ] * (Y[i, j] - p[i, j] - e[i, j])
    G[i, (ii + 1L):np] <- as.vector(glam)
  }
  Ainv <- solve(-object$hessian)                            # (-H)^{-1}
  Ainv %*% crossprod(G) %*% Ainv
}

#' Marginal Effects for a Mixed GCE Fit
#'
#' Partial derivatives of the signal probabilities with respect to each design
#' slot, holding \eqn{\lambda} and \eqn{\rho} fixed:
#' \deqn{\partial p_{ij}/\partial x_{ijk}
#'   = (\lambda_{kj}/\nu)\,\mathrm{Var}_\theta(s)_{ij},\qquad
#'   \mathrm{Var}_\theta(s)_{ij} = \sum_m s_m^2\theta_{ijm} - p_{ij}^2.}
#'
#' @details
#' With \code{se = TRUE} (and \code{average = TRUE}) standard errors for the
#' average marginal effects are attached. These are a \emph{sampling} quantity,
#' \strong{unrelated} to the Fano \code{\link{fano_bounds}} (which bound
#' classification error). \code{se_method = "sandwich"} (default) uses the robust
#' covariance \eqn{(-H)^{-1}\hat V(-H)^{-1}} of the dual multipliers, with
#' \eqn{\hat V=\sum_i g_i g_i'} the per-observation score outer products,
#' sandwiched with a numerical Jacobian of the average effects. This matches a
#' row bootstrap; the naive Hessian-inverse would overstate by ~7-9x for this
#' model (the support reparameterization plus per-observation \eqn{\rho} block
#' make \code{solve(-H)} a poor sampling-covariance estimate), so it is not used.
#' \code{se_method = "bootstrap"} resamples observations, refits, and takes the
#' across-resample SD.
#'
#' @param object A \code{mixed_gce} object.
#' @param average Logical (default \code{TRUE}). If \code{TRUE}, return the K x J
#'   matrix of marginal effects averaged over observations; if \code{FALSE}, the
#'   full N x J x K array of per-observation derivatives.
#' @param se Logical (default \code{FALSE}). If \code{TRUE} (requires
#'   \code{average = TRUE}), return a \code{margins_gce} object with estimates,
#'   standard errors, z-values, and p-values.
#' @param se_method \code{"sandwich"} (default) or \code{"bootstrap"}; see Details.
#' @param B Integer number of bootstrap resamples when
#'   \code{se_method = "bootstrap"} (default 500).
#' @param ... Unused.
#' @return With \code{se = FALSE}: a K x J matrix (average effects) or an
#'   N x J x K array. With \code{se = TRUE}: a \code{margins_gce} object.
#' @seealso \code{\link{fano_bounds}} for information-theoretic error bounds.
#' @export
margins.mixed_gce <- function(object, average = TRUE, se = FALSE,
                              se_method = c("sandwich", "bootstrap"),
                              B = 500L, ...) {
  theta <- object$theta; p <- object$p; lambda <- object$lambda
  s <- object$s; nu <- object$nu
  N <- object$N; J <- object$J; K <- object$K
  Es2 <- apply(sweep(theta, 3L, s^2, "*"), c(1L, 2L), sum)   # N x J: E_theta[s^2]
  vartheta <- Es2 - p^2                                      # N x J

  me_average <- function() {
    vbar <- colMeans(vartheta)
    out  <- matrix(0, K, J, dimnames = dimnames(lambda))
    for (j in seq_len(J)) out[, j] <- (lambda[, j] / nu) * vbar[j]
    out
  }

  if (!se) {
    if (average) return(me_average())
    ME <- array(0, dim = c(N, J, K),
                dimnames = list(NULL, colnames(lambda), rownames(lambda)))
    for (k in seq_len(K))
      for (j in seq_len(J)) ME[, j, k] <- (lambda[k, j] / nu) * vartheta[, j]
    return(ME)
  }

  # ---- se = TRUE: average marginal effects with standard errors -------------
  if (!average)
    stop("Standard errors are implemented only for the average effects ",
         "(average = TRUE).")
  se_method <- match.arg(se_method)
  ME <- me_average()
  theta0 <- object$theta0

  if (se_method == "sandwich") {
    par_hat <- c(object$rho, as.vector(lambda))
    np <- length(par_hat); h <- 1e-5
    Jac <- matrix(0, K * J, np)
    for (m in seq_len(np)) {
      e <- numeric(np); e[m] <- h
      Jac[, m] <- (.mixed_gce_me_vec(par_hat + e, object$X, s, theta0, nu, N, J, K) -
                   .mixed_gce_me_vec(par_hat - e, object$X, s, theta0, nu, N, J, K)) / (2 * h)
    }
    Vme <- Jac %*% .mixed_gce_sandwich(object) %*% t(Jac)
    SE  <- matrix(sqrt(pmax(diag(Vme), 0)), K, J, dimnames = dimnames(lambda))
  } else {
    MEb <- array(NA_real_, dim = c(K, J, B))
    for (b in seq_len(B)) {
      idx <- sample(N, N, replace = TRUE)
      fb  <- tryCatch(
        mixed_gce(object$y_mat[idx, , drop = FALSE], object$X[idx, , , drop = FALSE],
                  theta0 = theta0[idx, , , drop = FALSE],
                  w0 = object$w0[idx, , , drop = FALSE],
                  s = s, u = object$u, nu = nu),
        error = function(e) NULL)
      if (is.null(fb) || !fb$converged) next
      MEb[, , b] <- matrix(
        .mixed_gce_me_vec(c(fb$rho, as.vector(fb$lambda)),
                          object$X[idx, , , drop = FALSE], s,
                          theta0[idx, , , drop = FALSE], nu, N, J, K), K, J)
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
         se_method = se_method, n_boot = attr(SE, "n_boot")),
    class = "margins_gce"
  )
}

#' Fano Error Bounds for a Mixed GCE Fit
#'
#' For each observation (a row of \code{p_hat}, a distribution over the \eqn{J}
#' categories) the modal classifier predicts \eqn{\arg\max_j p_{ij}} with error
#' \eqn{pe_i = 1 - \max_j p_{ij}}, lower-bounded (Golan 2008, sec 7.5) by
#' \eqn{S(p_i) - \log 2/\log J}. This is an information-theoretic error bound,
#' \emph{not} a sampling standard error (for those see \code{\link{margins}} with
#' \code{se = TRUE}).
#'
#' @param object A \code{mixed_gce} object.
#' @param ... Unused.
#' @return A data frame with one row per observation (\code{p_max}, \code{pe},
#'   \code{H}, \code{S}, \code{pe_lower}) and an \code{"overall"} attribute.
#' @seealso \code{\link{fano_bounds}}, \code{\link{margins}}
#' @export
fano_bounds.mixed_gce <- function(object, ...) {
  .fano_row_bounds(object$p_hat, object$J)
}

# ===========================================================================
#  S3 methods
# ===========================================================================

#' @export
coef.mixed_gce <- function(object, ...) object$lambda

#' @export
fitted.mixed_gce <- function(object, ...) object$p_hat

#' @export
residuals.mixed_gce <- function(object, ...) object$y_mat - object$p_hat - object$e_hat

#' @export
vcov.mixed_gce <- function(object, ...) object$vcov

#' @export
print.mixed_gce <- function(x, digits = 4L, ...) {
  uniform <- max(abs(x$theta0 - 1 / x$M)) < 1e-10 && max(abs(x$w0 - 1 / x$H)) < 1e-10
  est <- if (uniform) "GME" else "GCE"
  cat(sprintf("\nMixed %s (doubly-reparameterized, dual)\n", est))
  cat(strrep("-", 52L), "\n")
  cat(sprintf("  N = %d   J = %d   K = %d   M = %d   H = %d   nu = %.3g\n",
              x$N, x$J, x$K, x$M, x$H, x$nu))
  cat(sprintf("  s support: [%.3g, %.3g]   u support: [%.3g, %.3g]\n",
              min(x$s), max(x$s), min(x$u), max(x$u)))
  cat(sprintf("  S(p)=%.4f  S(w)=%.4f  max|1-rowSums(p)|=%.2e  convergence=%s\n",
              x$S_p, x$S_w, max(abs(1 - rowSums(x$p))),
              if (x$converged) "yes" else "no"))
  cat("\nlambda (K x J):\n")
  print(round(x$lambda, digits))
  invisible(x)
}

#' @export
summary.mixed_gce <- function(object, digits = 4L, ...) {
  cat("\nMixed GCE (doubly-reparameterized, dual)\n")
  cat("Call: "); print(object$call); cat("\n")
  cat(sprintf("  N=%d  J=%d  K=%d  M=%d  H=%d  nu=%.3g\n",
              object$N, object$J, object$K, object$M, object$H, object$nu))
  cat(sprintf("  S(p)=%.4f  S(w)=%.4f  convergence=%s\n",
              object$S_p, object$S_w, if (object$converged) "yes" else "no"))
  fb <- attr(.fano_row_bounds(object$p_hat, object$J), "overall")
  cat(sprintf("  Fano (sec 7.5): mean modal error pe=%.4f  >= mean bound %.4f  [S(P)=%.4f]\n",
              fb["mean_pe"], fb["mean_pe_lower"], fb["S_system"]))
  for (j in seq_len(object$J)) {
    tab <- data.frame(
      lambda = object$lambda[, j],
      se     = object$se_lambda[, j],
      z      = ifelse(object$se_lambda[, j] > 0,
                      object$lambda[, j] / object$se_lambda[, j], NA_real_)
    )
    cat(sprintf("\n  Category '%s':\n", colnames(object$lambda)[j]))
    print(round(tab, digits))
  }
  cat(sprintf("\nAdding-up multipliers rho: range [%.3g, %.3g]\n",
              min(object$rho), max(object$rho)))
  cat("\nNote: coefficient SEs above are the (conservative) Hessian-inverse; for\n",
      "accurate marginal-effect SEs use margins(fit, se = TRUE). See ?margins.\n",
      sep = "")
  invisible(object)
}
