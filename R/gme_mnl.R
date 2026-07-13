# ============================================================
# gme_mnl.R — Multinomial IT estimators (ME and GME/GCE)
#             for unordered categorical response.
#
# Following Golan, Judge & Perloff (1996), JASA 91(434):841-853
# and the unified framework of Golan (2008), Chapter 6.
#
# PARAMETER CONVENTION
#   Lambda : K x J Lagrange multiplier matrix (col 1 = 0 by normalisation)
#   beta_j = -lambda_j   (ML-logit convention, eq. 28 of GJP 1996)
#   SE(beta) = SE(lambda)  (Jacobian of the sign flip is -I => identical)
#
# ME PRIMAL SOLUTIONS (uniform P0):
#   p_{ij}  = exp(-x_i' lambda_j) / Omega_i(lambda)
#
# GME/GCE PRIMAL SOLUTIONS (eqs. 25-26 of GJP 1996):
#   p_{ij}  = p0_{ij} * exp(-x_i' lambda_j)        / Omega_i(lambda)
#   w_{ijm} = w0_{ijm} * exp(-x_i' lambda_j * v_m)  / Psi_{ij}(lambda)
#   e_{ij}  = sum_m v_m * w_{ijm}
#
# DUAL OBJECTIVES (MINIMISED; positive-definite Hessian => convex):
#   ME:  L_ME(lam)  = sum Y*eta + sum_i log Omega_i
#   GME: L_GME(lam) = sum Y*eta + sum_i log Omega_i + sum_{ij} log Psi_{ij}
#
# GRADIENT (universal moment condition):
#   ME:  dL/d lambda_j = X'(Y - P)[, j]
#   GME: dL/d lambda_j = X'(Y - P - E)[, j]
#   (col 1 normalised to 0; free gradient for j = 2..J)
# ============================================================


# ---- Prior builders ---------------------------------------------------------

#' @keywords internal
.mnl_make_P0 <- function(p0, N, J) {
  if (is.null(p0))
    matrix(1 / J, nrow = N, ncol = J)
  else if (is.vector(p0) && length(p0) == J)
    matrix(p0, nrow = N, ncol = J, byrow = TRUE)
  else if (is.matrix(p0) && all(dim(p0) == c(N, J)))
    p0
  else
    stop("p0 must be NULL, a length-J vector, or an N x J matrix.")
}

#' @keywords internal
.mnl_make_W0 <- function(w0, N, J, M) {
  if (is.null(w0))
    array(1 / M, dim = c(N, J, M))
  else if (is.vector(w0) && length(w0) == M) {
    a <- array(0.0, c(N, J, M))
    for (m in seq_len(M)) a[, , m] <- w0[m]
    a
  } else if (is.array(w0) && all(dim(w0) == c(N, J, M)))
    w0
  else
    stop("w0 must be NULL, a length-M vector, or an N x J x M array.")
}

#' @keywords internal
.mnl_default_support <- function(N, M = 3L) {
  scale <- 1 / sqrt(N)
  seq(-scale, scale, length.out = M)
}


# ---- Vectorised probability computations ------------------------------------

# Compute N x J matrix P from Lambda (K x J) and X (N x K) with log prior.
# p_{ij} = p0_{ij} * exp(-eta_{ij}) / Omega_i  (log-sum-exp stabilised).

#' @keywords internal
.mnl_compute_P <- function(Lambda, X, logP0) {
  eta   <- X %*% Lambda                          # N x J
  arg_p <- -eta + logP0                          # log numerators (up to row const)
  m_p   <- apply(arg_p, 1, max)                  # row maxima for stability
  A     <- sweep(arg_p, 1, m_p, "-")
  expA  <- exp(A)
  expA / rowSums(expA)                           # N x J
}

# Compute P (N x J), W (N x J x M), and E (N x J) jointly.

#' @keywords internal
.mnl_compute_p_and_e <- function(Lambda, X, v, logP0, logW0) {
  N <- nrow(X); J <- ncol(Lambda); M <- length(v)
  eta <- X %*% Lambda                            # N x J

  # Signal probabilities (vectorised over observations)
  arg_p <- -eta + logP0
  m_p   <- apply(arg_p, 1, max)
  expA  <- exp(sweep(arg_p, 1, m_p, "-"))
  P     <- expA / rowSums(expA)

  # Noise weights and expected errors (loop over categories)
  W <- array(0.0, dim = c(N, J, M))
  E <- matrix(0.0, N, J)
  for (j in seq_len(J)) {
    Aj       <- -outer(eta[, j], v) + logW0[, j, ]   # N x M
    mj       <- apply(Aj, 1, max)
    eAj      <- exp(sweep(Aj, 1, mj, "-"))
    W[, j, ] <- eAj / rowSums(eAj)
    E[, j]   <- W[, j, ] %*% v
  }

  list(P = P, W = W, E = E)
}


# ---- Response extraction helper ---------------------------------------------

#' @keywords internal
.mnl_extract_response <- function(mf) {
  y_raw <- model.response(mf)
  if (is.factor(y_raw)) {
    levs  <- levels(y_raw)
    y_int <- as.integer(y_raw)
  } else if (is.character(y_raw)) {
    y_fac <- factor(y_raw)
    levs  <- levels(y_fac)
    y_int <- as.integer(y_fac)
  } else {
    y_int <- as.integer(y_raw)
    uvals <- sort(unique(y_int))
    levs  <- as.character(uvals)
    y_int <- match(y_int, uvals)
  }
  J <- length(levs)
  if (J < 2L) stop("Response must have at least 2 categories.")
  N <- length(y_int)
  Y <- matrix(0.0, N, J)
  Y[cbind(seq_len(N), y_int)] <- 1.0
  list(Y = Y, J = J, levs = levs)
}


# ---- ME dual objective and gradient -----------------------------------------

#' @keywords internal
.me_mnl_obj <- function(lambda_free, Y, X, logP0) {
  K <- ncol(X); J <- ncol(Y)
  Lambda        <- matrix(0.0, K, J)
  Lambda[, 2:J] <- matrix(lambda_free, K, J - 1)
  eta   <- X %*% Lambda
  arg_p <- -eta + logP0
  m_p   <- apply(arg_p, 1, max)
  sum(Y * eta) +
    sum(m_p + log(rowSums(exp(sweep(arg_p, 1, m_p, "-")))))
}

#' @keywords internal
.me_mnl_grad <- function(lambda_free, Y, X, logP0) {
  K <- ncol(X); J <- ncol(Y)
  Lambda        <- matrix(0.0, K, J)
  Lambda[, 2:J] <- matrix(lambda_free, K, J - 1)
  P <- .mnl_compute_P(Lambda, X, logP0)
  as.numeric(crossprod(X, (Y - P)[, 2:J, drop = FALSE]))
}


# ---- GME/GCE dual objective and gradient ------------------------------------

#' @keywords internal
.gme_mnl_obj <- function(lambda_free, Y, X, v, logP0, logW0) {
  K <- ncol(X); J <- ncol(Y)
  Lambda        <- matrix(0.0, K, J)
  Lambda[, 2:J] <- matrix(lambda_free, K, J - 1)
  eta <- X %*% Lambda

  # data term
  data_term <- sum(Y * eta)

  # signal partition  (log-sum-exp stabilised)
  arg_p       <- -eta + logP0
  m_p         <- apply(arg_p, 1, max)
  signal_term <- sum(m_p + log(rowSums(exp(sweep(arg_p, 1, m_p, "-")))))

  # noise partition
  noise_term <- 0.0
  for (j in seq_len(J)) {
    Aj         <- -outer(eta[, j], v) + logW0[, j, ]   # N x M
    mj         <- apply(Aj, 1, max)
    noise_term <- noise_term +
      sum(mj + log(rowSums(exp(sweep(Aj, 1, mj, "-")))))
  }

  data_term + signal_term + noise_term
}

#' @keywords internal
.gme_mnl_grad <- function(lambda_free, Y, X, v, logP0, logW0) {
  K <- ncol(X); J <- ncol(Y)
  Lambda        <- matrix(0.0, K, J)
  Lambda[, 2:J] <- matrix(lambda_free, K, J - 1)
  pe <- .mnl_compute_p_and_e(Lambda, X, v, logP0, logW0)
  as.numeric(crossprod(X, (Y - pe$P - pe$E)[, 2:J, drop = FALSE]))
}


# ---- Internal raw ME fitter (matrix interface; used for warm start) ---------

#' @keywords internal
.me_mnl_raw <- function(Y, X, init = NULL, control = list()) {
  N <- nrow(X); K <- ncol(X); J <- ncol(Y)
  logP0 <- matrix(log(1 / J), N, J)
  con <- list(maxit = 500L, reltol = 1e-10)
  con[names(control)] <- control
  if (is.null(init)) init <- rep(0.0, K * (J - 1))

  fit <- stats::optim(
    par     = init,
    fn      = .me_mnl_obj,
    gr      = .me_mnl_grad,
    Y = Y, X = X, logP0 = logP0,
    method  = "BFGS",
    control = list(maxit = con$maxit, reltol = con$reltol),
    hessian = TRUE
  )

  Lambda        <- matrix(0.0, K, J)
  Lambda[, 2:J] <- matrix(fit$par, K, J - 1)

  list(
    Lambda = Lambda,
    Beta   = -Lambda,
    P      = .mnl_compute_P(Lambda, X, logP0),
    fit    = fit,
    N = N, K = K, J = J
  )
}


# ===========================================================================
#  me_mnl()   —   ME estimator for multinomial choice (= multinomial logit)
# ===========================================================================

#' Maximum Entropy Multinomial Logit Estimator
#'
#' @description
#' Fits an unordered multinomial response model by Maximum Entropy (ME),
#' following the dual concentrated formulation of Golan, Judge & Perloff
#' (1996). With uniform prior \eqn{p_0 = 1/J} the ME solution coincides
#' with the standard multinomial logit maximum-likelihood estimator.
#'
#' @details
#' The dual objective (minimised, convex) is:
#' \deqn{
#'   L_{\mathrm{ME}}(\lambda) =
#'     \sum_{i,j} Y_{ij}\,\eta_{ij}
#'     + \sum_i \log\Omega_i(\lambda), \quad
#'   \eta_{ij} = x_i'\lambda_j, \quad
#'   \Omega_i = \sum_j \exp(-\eta_{ij}).
#' }
#' The normalisation \eqn{\lambda_{j=1}=0} is imposed, so \eqn{K(J-1)}
#' parameters are estimated. Coefficients are reported as
#' \eqn{\hat\beta_j = -\hat\lambda_j} to match the ML-logit convention
#' (eq. 28 of Golan, Judge & Perloff, 1996). Both the Hessian from
#' \code{optim} and the log-sum-exp trick are used for numerical stability.
#'
#' @param formula A symbolic formula. The left-hand side must be a factor
#'   (or coercible-to-integer) response; the right-hand side gives the
#'   regressors.
#' @param data A data frame containing the model variables.
#' @param init Optional numeric vector of length \eqn{K(J-1)} giving
#'   starting values for the free Lagrange multipliers. Defaults to zeros.
#' @param control Named list of solver settings: \code{maxit} (default
#'   \code{500}) and \code{reltol} (default \code{1e-10}).
#'
#' @return An object of class \code{c("infometrics_me_mnl", "infometrics")}
#'   with components:
#'   \describe{
#'     \item{\code{lambda}}{K x J matrix of Lagrange multipliers (col 1 = 0).}
#'     \item{\code{beta}}{K x J matrix, \code{-lambda}.}
#'     \item{\code{se_beta}}{K x J matrix of standard errors (from
#'       the inverse Hessian at the optimum).}
#'     \item{\code{vcov}}{Variance-covariance matrix of the free parameters,
#'       \eqn{H^{-1}}.}
#'     \item{\code{p_hat}}{N x J predicted probability matrix.}
#'     \item{\code{fitted.values}}{Same as \code{p_hat}.}
#'     \item{\code{loglik}}{Dual objective value at convergence.}
#'     \item{\code{S_p}}{Scalar normalised entropy of \code{p_hat} in
#'       \eqn{[0, 1]}.}
#'     \item{\code{misses}}{Number of misclassifications.}
#'     \item{\code{converged}}{Logical.}
#'     \item{\code{N, J, K}}{Dimensions.}
#'   }
#'
#' @references
#' Golan, A., Judge, G. and Perloff, J.M. (1996). A maximum entropy approach
#' to recovering information from multinomial response data.
#' \emph{Journal of the American Statistical Association}, \strong{91}(434),
#' 841-853.
#'
#' Golan, A. (2008). Information and Entropy Econometrics.
#' \emph{Foundations and Trends in Econometrics}, \strong{2}(1-2), 1-145.
#'
#' @seealso \code{\link{gme_mnl}} for the noise-extended GME/GCE estimator.
#'
#' @examples
#' set.seed(42)
#' n  <- 80
#' x  <- rnorm(n)
#' pr <- exp(cbind(0, 1 + 2 * x, -1 + x))
#' pr <- pr / rowSums(pr)
#' y  <- factor(apply(pr, 1, function(p) sample(3, 1, prob = p)))
#' df <- data.frame(y = y, x = x)
#' fit <- me_mnl(y ~ x, data = df)
#' print(fit)
#' coef(fit)
#'
#' @export
me_mnl <- function(formula, data, init = NULL, control = list()) {

  mc  <- match.call()
  mf  <- model.frame(formula, data)
  mt  <- attr(mf, "terms")
  X   <- model.matrix(mt, mf)
  N   <- nrow(X); K <- ncol(X)

  rsp  <- .mnl_extract_response(mf)
  Y    <- rsp$Y; J <- rsp$J; levs <- rsp$levs

  con <- list(maxit = 500L, reltol = 1e-10)
  con[names(control)] <- control

  if (!is.null(init) && length(init) != K * (J - 1))
    stop("init must have length K*(J-1) = ", K * (J - 1),
         ". Got length ", length(init), ".")

  # --- solve ---------------------------------------------------------------
  raw    <- .me_mnl_raw(Y, X, init = init, control = con)
  Lambda <- raw$Lambda; Beta <- raw$Beta; P <- raw$P; fit <- raw$fit

  # --- SE and vcov ---------------------------------------------------------
  vcov_free <- tryCatch(
    solve(fit$hessian),
    error = function(e) {
      warning("Hessian singular; SEs set to NA.")
      matrix(NA_real_, nrow(fit$hessian), ncol(fit$hessian))
    }
  )
  se_free <- sqrt(pmax(diag(vcov_free), 0))
  SE        <- matrix(0.0, K, J)
  SE[, 2:J] <- matrix(se_free, K, J - 1)

  # --- names ---------------------------------------------------------------
  rnames <- colnames(X)
  rownames(Lambda) <- rownames(Beta) <- rownames(SE) <- rnames
  colnames(Lambda) <- colnames(Beta) <- colnames(SE) <- levs
  colnames(P) <- levs

  # --- entropy and misclassification ---------------------------------------
  S_p    <- -sum(ifelse(P > 0, P * log(P), 0.0)) / (log(J) * N)
  misses <- sum(max.col(P) != max.col(Y))

  # per-observation signal entropy (N-vector; CLAUDE.md required field)
  H_signal <- -rowSums(ifelse(P > 0, P * log(P), 0.0))

  structure(
    list(
      call          = mc,
      terms         = mt,
      lambda        = Lambda,
      lambda_hat    = Lambda,   # CLAUDE.md alias
      beta          = Beta,
      se_lambda     = SE,
      se_beta       = SE,
      vcov          = vcov_free,
      p_hat         = P,
      fitted.values = P,
      y_mat         = Y,
      H_signal      = H_signal,
      S_p           = S_p,
      S             = S_p,      # CLAUDE.md alias
      loglik        = fit$value,
      objective     = fit$value,
      misses        = misses,
      converged     = (fit$convergence == 0L),
      iterations    = fit$counts,
      method        = "dual",
      N = N, J = J, K = K,
      levels        = levs,
      model         = mf,
      xlevels       = .getXlevels(mt, mf),
      contrasts     = attr(X, "contrasts")
    ),
    class = c("infometrics_me_mnl", "infometrics")
  )
}


# ===========================================================================
#  gme_mnl()   —   GME/GCE estimator for multinomial choice
# ===========================================================================

#' Generalized Maximum Entropy / Cross-Entropy Multinomial Estimator
#'
#' @description
#' Estimates an unordered multinomial response model via the Generalized
#' Maximum Entropy (GME) or Generalized Cross-Entropy (GCE) estimator of
#' Golan, Judge & Perloff (1996).
#'
#' With uniform priors (\code{p0 = NULL}, \code{w0 = NULL}) the estimator
#' maximises \eqn{H(p) + H(w)} (GME). With user-supplied priors it minimises
#' \eqn{D(p \| p_0) + D(w \| w_0)} (GCE). This mirrors \code{\link{gme}}:
#' uniform priors give GME, non-uniform priors give GCE — one function, not
#' two.
#'
#' @details
#' The dual objective (minimised, convex; Golan 2008, Eq. 6.5 generalised) is:
#' \deqn{
#'   L_{\mathrm{GME}}(\lambda) =
#'     \sum_{i,j} Y_{ij}\,\eta_{ij}
#'     + \sum_i \log\Omega_i(\lambda)
#'     + \sum_{i,j} \log\Psi_{ij}(\lambda),
#' }
#' where
#' \deqn{
#'   \Omega_i(\lambda) = \sum_j p_{0,ij}\exp(-\eta_{ij}), \quad
#'   \Psi_{ij}(\lambda) = \sum_m w_{0,ijm}\exp(-\eta_{ij}v_m).
#' }
#' Setting \eqn{p_0 = 1/J} and \eqn{w_0 = 1/M} recovers standard GME.
#' The Hessian is positive definite and BFGS converges to the unique
#' global minimum. An optional ME warm start (via \code{\link{me_mnl}})
#' accelerates convergence.
#'
#' @param formula A symbolic formula. The left-hand side must be a factor
#'   response.
#' @param data A data frame.
#' @param M Integer. Number of error support points. Ignored if \code{v}
#'   is supplied. Default \code{3}.
#' @param v Optional numeric vector of length M giving the error support
#'   (shared across all categories). If \code{NULL}, a symmetric M-point
#'   grid over \eqn{[-1/\sqrt{N},\, 1/\sqrt{N}]} is used (the data-scaled
#'   default of Golan, Judge & Perloff 1996).
#' @param p0 Signal prior: \code{NULL} (uniform, giving GME), a length-J
#'   probability vector, or an N x J matrix. Must be strictly positive
#'   with rows summing to 1. Supplying \code{p0} activates GCE.
#' @param w0 Error prior: \code{NULL} (uniform, giving GME), a length-M
#'   probability vector, or an N x J x M array. Must be strictly positive
#'   with \code{apply(w0, c(1,2), sum) == 1}. Supplying \code{w0}
#'   activates GCE.
#' @param init Optional numeric vector of length \eqn{K(J-1)} for starting
#'   values. If \code{NULL} and \code{warm_start = TRUE} an ME warm start
#'   is computed automatically.
#' @param warm_start Logical (default \code{TRUE}). Initialise from the ME
#'   (multinomial logit) solution for faster convergence.
#' @param control Named list: \code{maxit} (default \code{500}) and
#'   \code{reltol} (default \code{1e-10}).
#'
#' @return An object of class
#'   \code{c("infometrics_gme_mnl", "infometrics")} with components:
#'   \describe{
#'     \item{\code{lambda}}{K x J Lagrange multiplier matrix (col 1 = 0).}
#'     \item{\code{beta}}{K x J matrix, \code{-lambda}.}
#'     \item{\code{se_beta}}{K x J standard errors.}
#'     \item{\code{vcov}}{Variance-covariance of the \eqn{K(J-1)} free
#'       parameters, \eqn{H^{-1}}.}
#'     \item{\code{p_hat}}{N x J signal probability matrix.}
#'     \item{\code{w_hat}}{N x J x M noise weight array.}
#'     \item{\code{e}}{N x J recovered noise matrix.}
#'     \item{\code{S_p}}{Scalar normalised signal entropy in \eqn{[0,1]}.}
#'     \item{\code{S_w}}{Scalar normalised error entropy in \eqn{[0,1]}.}
#'     \item{\code{S_total}}{\code{S_p + S_w}.}
#'     \item{\code{KL_p}}{KL divergence \eqn{D(p \| p_0)} (0 for GME).}
#'     \item{\code{KL_w}}{KL divergence \eqn{D(w \| w_0)} (0 for GME).}
#'     \item{\code{misses}}{Number of misclassified observations.}
#'     \item{\code{converged}}{Logical.}
#'     \item{\code{N, J, K, M}}{Dimensions.}
#'   }
#'
#' @references
#' Golan, A., Judge, G. and Perloff, J.M. (1996). A maximum entropy approach
#' to recovering information from multinomial response data.
#' \emph{Journal of the American Statistical Association}, \strong{91}(434),
#' 841-853.
#'
#' Golan, A. (2008). Information and Entropy Econometrics.
#' \emph{Foundations and Trends in Econometrics}, \strong{2}(1-2), 1-145.
#'
#' @seealso \code{\link{me_mnl}} for the ME baseline (= multinomial logit);
#'   \code{\link{gme}} for the GME/GCE regression estimator.
#'
#' @examples
#' set.seed(42)
#' n  <- 80
#' x  <- rnorm(n)
#' pr <- exp(cbind(0, 1 + 2 * x, -1 + x))
#' pr <- pr / rowSums(pr)
#' y  <- factor(apply(pr, 1, function(p) sample(3, 1, prob = p)))
#' df <- data.frame(y = y, x = x)
#' fit <- gme_mnl(y ~ x, data = df)
#' print(fit)
#' coef(fit)
#'
#' @export
gme_mnl <- function(formula, data,
                    M = 3L, v = NULL,
                    p0 = NULL, w0 = NULL,
                    init = NULL, warm_start = TRUE,
                    control = list()) {

  mc  <- match.call()
  mf  <- model.frame(formula, data)
  mt  <- attr(mf, "terms")
  X   <- model.matrix(mt, mf)
  N   <- nrow(X); K <- ncol(X)

  rsp  <- .mnl_extract_response(mf)
  Y    <- rsp$Y; J <- rsp$J; levs <- rsp$levs

  con <- list(maxit = 500L, reltol = 1e-10)
  con[names(control)] <- control

  # ---- error support -------------------------------------------------------
  if (!is.null(v)) {
    v <- as.numeric(v); M <- length(v)
  } else {
    M <- as.integer(M)
    v <- .mnl_default_support(N, M)
  }
  if (M < 2L) stop("M must be >= 2.")

  # ---- priors --------------------------------------------------------------
  P0 <- .mnl_make_P0(p0, N, J)
  if (any(abs(rowSums(P0) - 1) > 1e-8))
    stop("Rows of p0 must sum to 1.")
  if (any(P0 <= 0))
    stop("p0 must be strictly positive.")

  W0 <- .mnl_make_W0(w0, N, J, M)
  if (any(abs(apply(W0, c(1, 2), sum) - 1) > 1e-8))
    stop("Each (i,j) fiber of w0 must sum to 1.")
  if (any(W0 <= 0))
    stop("w0 must be strictly positive.")

  logP0 <- log(P0); logW0 <- log(W0)

  # ---- initial values (ME warm start) -------------------------------------
  if (is.null(init)) {
    if (warm_start) {
      me_raw <- .me_mnl_raw(Y, X, control = list(maxit = 200L))
      # lambda_init = lambda_ME = -beta_ME
      init <- as.numeric(-me_raw$Beta[, 2:J, drop = FALSE])
    } else {
      init <- rep(0.0, K * (J - 1))
    }
  }
  if (length(init) != K * (J - 1))
    stop("init must have length K*(J-1) = ", K * (J - 1),
         ". Got length ", length(init), ".")

  # ---- optimise -----------------------------------------------------------
  fit <- stats::optim(
    par     = init,
    fn      = .gme_mnl_obj,
    gr      = .gme_mnl_grad,
    Y = Y, X = X, v = v, logP0 = logP0, logW0 = logW0,
    method  = "BFGS",
    control = list(maxit = con$maxit, reltol = con$reltol),
    hessian = TRUE
  )

  # ---- reconstruct primal solutions ---------------------------------------
  Lambda        <- matrix(0.0, K, J)
  Lambda[, 2:J] <- matrix(fit$par, K, J - 1)

  pe <- .mnl_compute_p_and_e(Lambda, X, v, logP0, logW0)
  P  <- pe$P; W <- pe$W; E <- pe$E
  Beta <- -Lambda

  # ---- SE and vcov --------------------------------------------------------
  vcov_free <- tryCatch(
    solve(fit$hessian),
    error = function(e) {
      warning("Hessian singular; SEs set to NA.")
      matrix(NA_real_, nrow(fit$hessian), ncol(fit$hessian))
    }
  )
  se_free   <- sqrt(pmax(diag(vcov_free), 0))
  SE        <- matrix(0.0, K, J)
  SE[, 2:J] <- matrix(se_free, K, J - 1)

  # ---- names --------------------------------------------------------------
  rnames <- colnames(X)
  rownames(Lambda) <- rownames(Beta) <- rownames(SE) <- rnames
  colnames(Lambda) <- colnames(Beta) <- colnames(SE) <- levs
  colnames(P) <- colnames(E) <- levs

  # ---- entropy measures ---------------------------------------------------
  Plog   <- ifelse(P > 0, P * log(P), 0.0)
  S_p    <- -sum(Plog) / (log(J) * N)

  Wf     <- as.numeric(W)
  Wlog   <- ifelse(Wf > 0, Wf * log(Wf), 0.0)
  S_w    <- -sum(Wlog) / (log(M) * N * J)

  KL_p   <- sum(ifelse(P > 0, P * (log(P) - logP0), 0.0))
  KL_w   <- sum(ifelse(Wf > 0, Wf * (log(Wf) - as.numeric(logW0)), 0.0))

  misses <- sum(max.col(P) != max.col(Y))

  # per-observation entropy (N-vector; CLAUDE.md required field)
  H_signal <- -rowSums(Plog)
  H_error  <- -rowSums(matrix(Wlog, N, J * M))[seq_len(N)]  # scalar total
  # simpler scalar total error entropy
  H_error_total <- -sum(Wlog)

  structure(
    list(
      call          = mc,
      terms         = mt,
      lambda        = Lambda,
      lambda_hat    = Lambda,      # CLAUDE.md alias
      beta          = Beta,
      se_lambda     = SE,
      se_beta       = SE,
      vcov          = vcov_free,
      p_hat         = P,
      w_hat         = W,
      fitted.values = P,
      e             = E,
      y_mat         = Y,
      P0            = P0,
      W0            = W0,
      v             = v,
      M             = M,
      H_signal      = H_signal,
      H_error       = H_error_total,
      loglik        = fit$value,
      objective     = fit$value,
      S_p           = S_p,
      S_w           = S_w,
      S_total       = S_p + S_w,
      S             = S_p + S_w,   # CLAUDE.md alias
      KL_p          = KL_p,
      KL_w          = KL_w,
      misses        = misses,
      hessian       = fit$hessian,
      converged     = (fit$convergence == 0L),
      iterations    = fit$counts,
      method        = "dual",
      N = N, J = J, K = K, M = M,
      levels        = levs,
      model         = mf,
      xlevels       = .getXlevels(mt, mf),
      contrasts     = attr(X, "contrasts")
    ),
    class = c("infometrics_gme_mnl", "infometrics")
  )
}


# ===========================================================================
#  S3 methods — infometrics_me_mnl
# ===========================================================================

#' @export
print.infometrics_me_mnl <- function(x, digits = 4L, ...) {
  cat("\nMaximum Entropy Multinomial Logit (ME-MNL)\n")
  cat("-------------------------------------------\n")
  cat(sprintf("  N = %d   J = %d   K = %d\n", x$N, x$J, x$K))
  cat(sprintf("  L_ME       : %.4f\n", x$loglik))
  cat(sprintf("  S(p_hat)   : %.4f\n", x$S_p))
  cat(sprintf("  misses     : %d / %d (%.1f%%)\n",
              x$misses, x$N, 100 * x$misses / x$N))
  cat(sprintf("  convergence: %s\n", if (x$converged) "yes" else "no"))
  cat("\nBeta = -Lambda (ML-logit convention; baseline = ",
      x$levels[1], "):\n", sep = "")
  print(round(x$beta, digits))
  invisible(x)
}

#' @export
summary.infometrics_me_mnl <- function(object, digits = 4L, ...) {
  J    <- object$J
  tabs <- lapply(seq_len(J)[-1], function(j) {
    data.frame(
      beta   = object$beta[, j],
      se     = object$se_beta[, j],
      lambda = object$lambda[, j],
      t_stat = ifelse(object$se_beta[, j] > 0,
                      object$beta[, j] / object$se_beta[, j],
                      NA_real_)
    )
  })
  names(tabs) <- object$levels[-1]

  cat("\nMaximum Entropy Multinomial Logit (ME-MNL)\n")
  cat("Call: "); print(object$call); cat("\n")
  cat(sprintf(
    "  N=%d  J=%d  K=%d  L=%.4f  S(p)=%.4f  misses=%d/%d\n",
    object$N, object$J, object$K,
    object$loglik, object$S_p, object$misses, object$N
  ))
  for (nm in names(tabs)) {
    cat(sprintf("\n  Category '%s' (vs. baseline '%s'):\n",
                nm, object$levels[1]))
    print(round(tabs[[nm]], digits))
  }
  invisible(object)
}

#' @export
coef.infometrics_me_mnl <- function(object, ...) {
  object$beta
}

#' @export
fitted.infometrics_me_mnl <- function(object, ...) {
  object$fitted.values
}

#' @export
residuals.infometrics_me_mnl <- function(object, ...) {
  object$y_mat - object$p_hat
}


# ===========================================================================
#  S3 methods — infometrics_gme_mnl
# ===========================================================================

#' @export
print.infometrics_gme_mnl <- function(x, digits = 4L, ...) {
  uniform_p  <- max(abs(x$P0 - 1 / x$J)) < 1e-10
  uniform_w  <- max(abs(x$W0 - 1 / x$M)) < 1e-10
  estimator  <- if (uniform_p && uniform_w) "GME" else "GCE"
  long_name  <- c(GME = "Generalized Maximum Entropy",
                  GCE = "Generalized Cross-Entropy")[estimator]

  cat(sprintf("\n%s Multinomial (%s)\n", long_name, estimator))
  cat(strrep("-", 50L), "\n")
  cat(sprintf("  N = %d   J = %d   K = %d   M = %d\n", x$N, x$J, x$K, x$M))
  cat(sprintf("  prior p0   : %s\n",
              if (uniform_p) "uniform (1/J)" else "non-uniform"))
  cat(sprintf("  prior w0   : %s\n",
              if (uniform_w) "uniform (1/M)" else "non-uniform"))
  cat(sprintf("  support v  : [%.4f, %.4f]  M = %d\n",
              min(x$v), max(x$v), x$M))
  cat(sprintf("  L_GME      : %.4f\n", x$loglik))
  cat(sprintf("  S(p_hat)   : %.4f   KL(p||p0) = %.4f\n", x$S_p, x$KL_p))
  cat(sprintf("  S(w_hat)   : %.4f   KL(w||w0) = %.4f\n", x$S_w, x$KL_w))
  cat(sprintf("  misses     : %d / %d (%.1f%%)\n",
              x$misses, x$N, 100 * x$misses / x$N))
  cat(sprintf("  convergence: %s\n", if (x$converged) "yes" else "no"))
  cat(sprintf("\nBeta = -Lambda (baseline = '%s'):\n", x$levels[1]))
  print(round(x$beta, digits))
  invisible(x)
}

#' @export
summary.infometrics_gme_mnl <- function(object, digits = 4L, ...) {
  J    <- object$J
  tabs <- lapply(seq_len(J)[-1], function(j) {
    data.frame(
      beta   = object$beta[, j],
      se     = object$se_beta[, j],
      lambda = object$lambda[, j],
      t_stat = ifelse(object$se_beta[, j] > 0,
                      object$beta[, j] / object$se_beta[, j],
                      NA_real_)
    )
  })
  names(tabs) <- object$levels[-1]

  uniform_p <- max(abs(object$P0 - 1 / object$J)) < 1e-10
  uniform_w <- max(abs(object$W0 - 1 / object$M)) < 1e-10
  estimator <- if (uniform_p && uniform_w) "GME" else "GCE"
  long_name <- c(GME = "Generalized Maximum Entropy",
                 GCE = "Generalized Cross-Entropy")[estimator]

  cat(sprintf("\n%s Multinomial (%s)\n", long_name, estimator))
  cat("Call: "); print(object$call); cat("\n")
  cat(sprintf("  N=%d  J=%d  K=%d  M=%d  L=%.4f\n",
              object$N, object$J, object$K, object$M, object$loglik))
  cat(sprintf(
    "  S(p)=%.4f  S(w)=%.4f  S_total=%.4f  misses=%d/%d\n",
    object$S_p, object$S_w, object$S_total,
    object$misses, object$N
  ))
  for (nm in names(tabs)) {
    cat(sprintf("\n  Category '%s' (vs. baseline '%s'):\n",
                nm, object$levels[1]))
    print(round(tabs[[nm]], digits))
  }
  invisible(object)
}

#' @export
coef.infometrics_gme_mnl <- function(object, ...) {
  object$beta
}

#' @export
fitted.infometrics_gme_mnl <- function(object, ...) {
  object$fitted.values
}

#' @export
residuals.infometrics_gme_mnl <- function(object, ...) {
  object$y_mat - object$p_hat - object$e
}
