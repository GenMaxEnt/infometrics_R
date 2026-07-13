# ============================================================
# gce_clogit.R — Conditional/Mixed Logit GCE estimator
#
# Implements the GCE extension of the conditional logit model
# following Golan, Judge & Perloff (1996) and the mixed logit
# specification of the mlogit package (Croissant, 2020).
#
# FORMULA INTERFACE (mlogit-style):
#   choice ~ choice_specific | individual_specific
#
#   Left of |  : variables that vary by BOTH individual AND alternative
#                -> one shared beta per variable ("generic coefficients")
#   Right of | : variables that vary by individual ONLY
#                -> (J-1) alternative-specific gammas per variable
#   Alternative-specific constants (J-1) are always added automatically.
#   Reference alternative = first level of the choice factor.
#
# DATA FORMAT: long format, one row per (individual, alternative) pair,
# sorted consistently by (individual-id, alternative). The formula LHS
# column contains TRUE/1 for the chosen row, FALSE/0 otherwise.
#
# PARAMETER CONVENTION:
#   lambda : K-vector of Lagrange multipliers (all free, no normalization
#            since alternative specificity is embedded in Z)
#   beta   = -lambda   (ML-logit / clogit convention; reported by default)
#   SE(beta) = SE(lambda)   [Jacobian of sign flip is -I => identical]
#
# UTILITY: V_ij = z_ij' lambda   where z_ij = Z[i, j, ] is a K-vector
#
# PRIMAL SOLUTIONS:
#   p_hat_ij  = p0_ij * exp(-V_ij) / Omega_i(lambda)
#   w_hat_ijm = w0_ijm * exp(-V_ij * v_m) / Psi_ij(lambda)
#   e_hat_ij  = sum_m v_m * w_hat_ijm
#
# DUAL OBJECTIVE (MINIMISED; convex; Golan 2008 Eq. 6.5 generalised):
#   L_GCE(lambda) = sum_{ij} y_ij V_ij
#                 + sum_i log Omega_i^0(lambda)
#                 + sum_{ij} log Psi_ij^0(lambda)
#
# GRADIENT (universal moment condition):
#   dL / d lambda = sum_{ij} (y_ij - p_hat_ij - e_hat_ij) z_ij
#
# Uniform priors (p0 = NULL, w0 = NULL) -> GME-clogit (max entropy).
# User-supplied priors -> GCE-clogit (min cross-entropy).
# ============================================================


# ---- Formula parser ---------------------------------------------------------

#' @keywords internal
.parse_clogit_formula <- function(formula) {
  fchar <- paste(deparse(formula, width.cutoff = 500L), collapse = " ")

  parts_lhs <- strsplit(fchar, "~", fixed = TRUE)[[1]]
  if (length(parts_lhs) != 2L)
    stop("Formula must be of form  choice ~ ... | ...")
  lhs <- trimws(parts_lhs[1])
  rhs <- trimws(parts_lhs[2])

  rhs_parts <- trimws(strsplit(rhs, "|", fixed = TRUE)[[1]])
  if (length(rhs_parts) > 3L)
    stop("Formula has more than 3 '|' sections (max 3).")
  rhs_parts[rhs_parts == ""] <- "0"

  list(
    response   = lhs,
    choice_sp  = as.formula(paste("~", rhs_parts[1])),
    indiv_sp   = if (length(rhs_parts) >= 2L)
                   as.formula(paste("~", rhs_parts[2]))
                 else as.formula("~ 0"),
    alt_sp     = if (length(rhs_parts) >= 3L)
                   as.formula(paste("~", rhs_parts[3]))
                 else as.formula("~ 0")
  )
}


# ---- Design matrix builder --------------------------------------------------
#
# Returns Z (N x J x K), Y (N x J), alt_levels, parameter names.
#
# K = (J-1) ASCs + K_cs generic coefs + (J-1)*K_is alt-specific coefs
#
# Z layout (column order in the K dimension):
#   [ alpha_2, ..., alpha_J,          <- J-1 ASCs
#     beta_1, ..., beta_{K_cs},       <- choice-specific generic coefs
#     gamma_{2,1}, ..., gamma_{J,1},  <- alt-specific coefs for var 1
#     ...                             <- (repeated for each indiv-spec var)
#   ]

#' @keywords internal
.build_clogit_design <- function(data, formula, alt, id) {

  parts      <- .parse_clogit_formula(formula)
  alt_levels <- unique(data[[alt]])
  ids        <- unique(data[[id]])
  J <- length(alt_levels)
  N <- length(ids)

  if (nrow(data) != N * J)
    stop(sprintf(
      "Data has %d rows but N x J = %d x %d = %d expected. ",
      nrow(data), N, J, N * J),
      "Ensure every individual has exactly one row per alternative.")

  # Sort: id primary, alt secondary (matching alt_levels order)
  data <- data[order(data[[id]], match(data[[alt]], alt_levels)), , drop = FALSE]

  # --- Choice-specific variables (generic coefs) ---
  X_cs <- if (length(all.vars(parts$choice_sp)) > 0)
    model.matrix(parts$choice_sp, data)[, -1L, drop = FALSE]
  else
    matrix(0.0, nrow(data), 0L)
  K_cs <- ncol(X_cs)

  # --- Individual-specific variables (alt-specific coefs) ---
  # One row per individual suffices (all rows of same id are identical)
  data_ind <- data[!duplicated(data[[id]]), , drop = FALSE]
  X_is <- if (length(all.vars(parts$indiv_sp)) > 0)
    model.matrix(parts$indiv_sp, data_ind)[, -1L, drop = FALSE]
  else
    matrix(0.0, N, 0L)
  K_is <- ncol(X_is)

  K <- (J - 1L) + K_cs + (J - 1L) * K_is
  Z <- array(0.0, c(N, J, K))
  cnames <- character(K)
  pos <- 1L

  # Alternative-specific constants (J-1; reference = alt_levels[1])
  for (j in 2:J) {
    Z[, j, pos] <- 1.0
    cnames[pos] <- paste0("(Intercept):", alt_levels[j])
    pos <- pos + 1L
  }

  # Choice-specific variables: reshape (N*J) rows -> N x J matrix per var
  if (K_cs > 0L) {
    for (k in seq_len(K_cs)) {
      vals <- matrix(X_cs[, k], nrow = N, ncol = J, byrow = TRUE)
      Z[, , pos] <- vals
      cnames[pos] <- colnames(X_cs)[k]
      pos <- pos + 1L
    }
  }

  # Individual-specific variables: place X_is[, k] in col j of Z for j >= 2
  if (K_is > 0L) {
    for (k in seq_len(K_is)) {
      vname <- colnames(X_is)[k]
      for (j in 2:J) {
        Z[, j, pos] <- X_is[, k]
        cnames[pos] <- paste0(vname, ":", alt_levels[j])
        pos <- pos + 1L
      }
    }
  }

  # Response: Y is N x J indicator (chosen = 1)
  choice_vec <- as.logical(data[[parts$response]])
  Y <- matrix(0.0, N, J)
  for (i in seq_len(N)) {
    block  <- choice_vec[(i - 1L) * J + seq_len(J)]
    chosen <- which(block)
    if (length(chosen) != 1L)
      stop("Individual ", ids[i], " has ", length(chosen),
           " chosen row(s); exactly 1 expected.")
    Y[i, chosen] <- 1.0
  }

  list(Z = Z, Y = Y, names = cnames, alt_levels = alt_levels,
       N = N, J = J, K = K)
}


# ---- Compute utilities V (N x J) from Z and lambda -------------------------

#' @keywords internal
.clogit_utilities <- function(lambda, Z) {
  J <- dim(Z)[2]
  # V[, j] = Z[, j, ] %*% lambda  for each j
  vapply(seq_len(J), function(j) drop(Z[, j, ] %*% lambda), numeric(dim(Z)[1]))
}


# ---- Compute P (N x J), W (N x J x M), E (N x J) from V -------------------

#' @keywords internal
.clogit_compute_pe <- function(V, v, logP0, logW0) {
  N <- nrow(V); J <- ncol(V); M <- length(v)

  # Signal probabilities (vectorised)
  arg_p <- -V + logP0
  m_p   <- apply(arg_p, 1L, max)
  expA  <- exp(sweep(arg_p, 1L, m_p, "-"))
  P     <- expA / rowSums(expA)

  # Noise weights and expected errors
  W <- array(0.0, dim = c(N, J, M))
  E <- matrix(0.0, N, J)
  for (j in seq_len(J)) {
    Aj       <- -outer(V[, j], v) + logW0[, j, ]   # N x M
    mj       <- apply(Aj, 1L, max)
    eAj      <- exp(sweep(Aj, 1L, mj, "-"))
    W[, j, ] <- eAj / rowSums(eAj)
    E[, j]   <- W[, j, ] %*% v
  }

  list(P = P, W = W, E = E)
}


# ---- Analytic Hessian of the GCE dual objective ----------------------------
#
#  H = Σ_i [ Z_i' diag(p_i) Z_i − (Z_i' p_i)(Z_i' p_i)' ]   (signal)
#    + Σ_{ij} Var_{w_ij}(v) · z_ij z_ij'                      (noise)
#
#  Derivation:
#    Signal: d² log Omega_i / d λ d λ' = Var_{p_i}[z_i]
#            = Z_i' diag(p_i) Z_i − (Z_i' p_i)(Z_i' p_i)'
#    Noise:  d² log Psi_ij / d λ d λ' = Var_{w_ij}[v] · z_ij z_ij'
#            where Var_{w_ij}[v] = Σ_m w_ijm v_m² − (Σ_m w_ijm v_m)²
#
#  Both terms are positive semidefinite; their sum is always PSD,
#  guaranteeing solve(H) produces valid (non-negative) diagonal entries.

#' @keywords internal
.clogit_analytic_hessian <- function(Z, P, W, v) {
  N  <- dim(Z)[1L]; J <- dim(Z)[2L]; K <- dim(Z)[3L]
  v2 <- v^2
  H  <- matrix(0.0, K, K)

  for (i in seq_len(N)) {
    Zi  <- Z[i, , ]                    # J × K
    pi  <- P[i, ]                      # J-vector
    Ztp <- drop(crossprod(Zi, pi))     # K-vector: Z_i' p_i

    # Signal contribution: Z_i' diag(p_i) Z_i − (Z_i' p_i)(Z_i' p_i)'
    H <- H + crossprod(Zi, pi * Zi) - outer(Ztp, Ztp)

    # Noise contribution: Σ_j Var_{w_ij}(v) · z_ij z_ij'
    for (j in seq_len(J)) {
      wij    <- W[i, j, ]
      eij    <- sum(wij * v)
      var_ij <- sum(wij * v2) - eij^2
      zij    <- Z[i, j, ]
      H      <- H + var_ij * outer(zij, zij)
    }
  }
  H
}


# ---- ME conditional logit (standard CL, no noise) — warm start -------------

#' @keywords internal
.clogit_me_obj <- function(lambda, Y, Z, logP0) {
  V     <- .clogit_utilities(lambda, Z)
  arg_p <- -V + logP0
  m_p   <- apply(arg_p, 1L, max)
  sum(Y * V) +
    sum(m_p + log(rowSums(exp(sweep(arg_p, 1L, m_p, "-")))))
}

#' @keywords internal
.clogit_me_grad <- function(lambda, Y, Z, logP0) {
  J     <- dim(Z)[2]
  V     <- .clogit_utilities(lambda, Z)
  arg_p <- -V + logP0
  m_p   <- apply(arg_p, 1L, max)
  expA  <- exp(sweep(arg_p, 1L, m_p, "-"))
  P     <- expA / rowSums(expA)
  res   <- Y - P
  Reduce("+", lapply(seq_len(J), function(j) colSums(res[, j] * Z[, j, ])))
}

#' @keywords internal
.clogit_me_raw <- function(Y, Z, init = NULL, control = list()) {
  N <- dim(Z)[1]; K <- dim(Z)[3]
  J <- dim(Z)[2]
  logP0 <- matrix(log(1 / J), N, J)
  con <- list(maxit = 500L, reltol = 1e-10); con[names(control)] <- control
  if (is.null(init)) init <- rep(0.0, K)

  fit <- stats::optim(
    par     = init,
    fn      = .clogit_me_obj,
    gr      = .clogit_me_grad,
    Y = Y, Z = Z, logP0 = logP0,
    method  = "BFGS",
    control = list(maxit = con$maxit, reltol = con$reltol)
  )
  list(lambda = fit$par, beta = -fit$par, fit = fit)
}


# ---- GCE-clogit dual objective and gradient ---------------------------------

#' @keywords internal
.clogit_obj <- function(lambda, Y, Z, v, logP0, logW0) {
  J <- dim(Z)[2]
  V <- .clogit_utilities(lambda, Z)   # N x J

  # Data term
  data_term <- sum(Y * V)

  # Signal term: sum_i log Omega_i (vectorised)
  arg_p       <- -V + logP0
  m_p         <- apply(arg_p, 1L, max)
  signal_term <- sum(m_p + log(rowSums(exp(sweep(arg_p, 1L, m_p, "-")))))

  # Noise term: sum_{ij} log Psi_ij (loop over j, vectorised over i)
  noise_term <- 0.0
  for (j in seq_len(J)) {
    Aj <- -outer(V[, j], v) + logW0[, j, ]   # N x M
    mj <- apply(Aj, 1L, max)
    noise_term <- noise_term +
      sum(mj + log(rowSums(exp(sweep(Aj, 1L, mj, "-")))))
  }

  data_term + signal_term + noise_term
}

#' @keywords internal
.clogit_grad <- function(lambda, Y, Z, v, logP0, logW0) {
  J <- dim(Z)[2]
  V  <- .clogit_utilities(lambda, Z)
  pe <- .clogit_compute_pe(V, v, logP0, logW0)

  res <- Y - pe$P - pe$E    # N x J residual matrix
  Reduce("+", lapply(seq_len(J), function(j) colSums(res[, j] * Z[, j, ])))
}


# ===========================================================================
#  gce_clogit()  —  Conditional / Mixed Logit GCE estimator
# ===========================================================================

#' Conditional / Mixed Logit GCE Estimator
#'
#' @description
#' Estimates a conditional logit or mixed logit model via the Generalized
#' Cross-Entropy (GCE) principle, following Golan, Judge & Perloff (1996)
#' and the mixed logit specification of Croissant (2020).
#'
#' The estimator unifies three variable types through a single
#' \code{mlogit}-style formula:
#' \itemize{
#'   \item \strong{Choice-specific variables} (left of \code{|}): attributes
#'     that vary by both individual and alternative (e.g. price, travel time).
#'     One shared coefficient per variable.
#'   \item \strong{Individual-specific variables} (right of \code{|}):
#'     attributes that vary only by individual (e.g. income). Separate
#'     coefficient per non-reference alternative.
#'   \item \strong{Alternative-specific constants}: added automatically
#'     for every alternative except the reference (first level).
#' }
#'
#' @details
#' The utility for individual \eqn{i} choosing alternative \eqn{j} is
#' \eqn{V_{ij} = z_{ij}'\lambda}, where \eqn{z_{ij}} is a \eqn{K}-vector
#' assembled from the three variable types above, and \eqn{\lambda} is a
#' single \eqn{K}-vector of Lagrange multipliers.  Coefficients are reported
#' as \eqn{\hat\beta = -\hat\lambda} (ML-clogit convention).
#'
#' The dual objective (minimised, convex; Golan 2008, Eq. 6.5 generalised):
#' \deqn{
#'   L_{\mathrm{GCE}}(\lambda) = \sum_{ij} y_{ij} V_{ij}
#'   + \sum_i \log\Omega_i(\lambda)
#'   + \sum_{ij} \log\Psi_{ij}(\lambda),
#' }
#' where \eqn{\Omega_i = \sum_j p_{0,ij}\exp(-V_{ij})} and
#' \eqn{\Psi_{ij} = \sum_m w_{0,ijm}\exp(-V_{ij}v_m)}.
#' Uniform priors (\code{p0 = NULL}, \code{w0 = NULL}) give the
#' maximum-entropy (GME) version; user-supplied priors give GCE.
#'
#' @param formula A two-part formula of the form
#'   \code{chosen ~ choice_vars | indiv_vars}.
#'   \code{chosen} is a logical/0-1 column in the long-format data that
#'   equals \code{TRUE} for the row corresponding to the selected
#'   alternative.  Either part may be \code{0} (no variables of that type).
#'   A third part \code{| alt_vars} is parsed but currently ignored.
#' @param data A data frame in \strong{long format}: one row per
#'   (individual, alternative) pair.  Must contain columns named by
#'   \code{alt} and \code{id}.
#' @param alt Character. Name of the column that identifies the alternative
#'   (e.g. \code{"mode"}).
#' @param id  Character. Name of the column that identifies the individual
#'   (e.g. \code{"chid"}).
#' @param M   Integer. Number of error support points (default \code{3}).
#'   Ignored if \code{v} is supplied.
#' @param v   Optional numeric vector of length \code{M} giving the error
#'   support.  If \code{NULL} a symmetric \code{M}-point grid over
#'   \eqn{[-1/\sqrt{N},\, 1/\sqrt{N}]} is used.
#' @param p0  Optional signal prior.  \code{NULL} (uniform, GME), a
#'   length-\eqn{J} probability vector, or an \eqn{N \times J} matrix.
#'   Must be strictly positive with rows summing to 1.
#' @param w0  Optional error prior.  \code{NULL} (uniform, GME), a
#'   length-\eqn{M} probability vector, or an
#'   \eqn{N \times J \times M} array.
#' @param init Optional numeric vector of length \eqn{K} for starting
#'   values.  If \code{NULL} and \code{warm_start = TRUE} a standard
#'   conditional logit (ME) solution is used.
#' @param warm_start Logical (default \code{TRUE}).  Initialise from the
#'   maximum-likelihood conditional logit solution for faster BFGS
#'   convergence.
#' @param control Named list: \code{maxit} (default \code{500}) and
#'   \code{reltol} (default \code{1e-10}).
#'
#' @return An object of class \code{c("infometrics_clogit", "infometrics")}
#'   with components:
#'   \describe{
#'     \item{\code{beta}}{Named \eqn{K}-vector of estimated coefficients,
#'       \eqn{\hat\beta = -\hat\lambda}.}
#'     \item{\code{lambda}}{Named \eqn{K}-vector of Lagrange multipliers.}
#'     \item{\code{se_beta}}{Named \eqn{K}-vector of standard errors.}
#'     \item{\code{vcov}}{\eqn{K \times K} variance-covariance matrix,
#'       \eqn{H^{-1}}.}
#'     \item{\code{p_hat}}{\eqn{N \times J} matrix of predicted choice
#'       probabilities.}
#'     \item{\code{w_hat}}{\eqn{N \times J \times M} noise weight array.}
#'     \item{\code{e}}{\eqn{N \times J} recovered noise matrix.}
#'     \item{\code{loglik}}{Log-likelihood at convergence
#'       (\eqn{-L_{\mathrm{GCE}}}).}
#'     \item{\code{S_p}}{Scalar normalised signal entropy in \eqn{[0,1]}.}
#'     \item{\code{S_w}}{Scalar normalised error entropy in \eqn{[0,1]}.}
#'     \item{\code{KL_p}}{KL divergence \eqn{D(p \| p_0)}.}
#'     \item{\code{KL_w}}{KL divergence \eqn{D(w \| w_0)}.}
#'     \item{\code{misses}}{Number of misclassified individuals.}
#'     \item{\code{converged}}{Logical.}
#'     \item{\code{N, J, K, M}}{Dimensions.}
#'     \item{\code{alt_levels}}{Character vector of alternative labels.}
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
#' Croissant, Y. (2020). Estimation of random utility models in R: The mlogit
#' package. \emph{Journal of Statistical Software}, \strong{95}(11), 1-41.
#'
#' @seealso \code{\link{gme_mnl}} for the GME/GCE multinomial estimator
#'   (wide-format, no choice-specific variables).
#'
#' @examples
#' ## Fishing mode choice data (from the mlogit package)
#' ## Recreated as a small subset for the example
#' set.seed(42)
#' n_ind <- 30L
#' alts  <- c("beach", "boat", "charter", "pier")
#' J     <- length(alts)
#' N     <- n_ind
#'
#' # Simulate long-format data
#' df_long <- data.frame(
#'   id    = rep(seq_len(N), each = J),
#'   alt   = rep(alts, times = N),
#'   price = round(runif(N * J, 50, 300)),
#'   catch = round(runif(N * J, 0.1, 1.5), 2),
#'   income = rep(round(runif(N, 20, 80)), each = J)
#' )
#' # Simulate choices via a simple logit
#' eta   <- with(df_long, -0.01 * price + 0.5 * catch)
#' prob  <- tapply(exp(eta), df_long$id, function(x) x / sum(x))
#' chosen <- unlist(lapply(seq_len(N), function(i)
#'   seq_len(J) == sample(J, 1, prob = prob[[i]])))
#' df_long$chosen <- chosen
#'
#' fit <- gce_clogit(chosen ~ price + catch | income,
#'                   data = df_long, alt = "alt", id = "id")
#' print(fit)
#' coef(fit)
#'
#' @export
gce_clogit <- function(formula, data,
                       alt, id,
                       M = 3L, v = NULL,
                       p0 = NULL, w0 = NULL,
                       init = NULL, warm_start = TRUE,
                       control = list()) {

  mc  <- match.call()

  # ---- build design --------------------------------------------------------
  d <- .build_clogit_design(data, formula, alt = alt, id = id)
  Z <- d$Z; Y <- d$Y; N <- d$N; J <- d$J; K <- d$K

  # ---- control defaults ----------------------------------------------------
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

  # ---- initial values (ME/CL warm start) -----------------------------------
  if (is.null(init)) {
    if (warm_start) {
      me_raw <- .clogit_me_raw(Y, Z, control = list(maxit = 200L))
      init   <- me_raw$lambda   # lambda_ME is a good start for GCE
    } else {
      init <- rep(0.0, K)
    }
  }
  if (length(init) != K)
    stop("init must have length K = ", K, ". Got length ", length(init), ".")

  # ---- optimise ------------------------------------------------------------
  fit <- stats::optim(
    par     = init,
    fn      = .clogit_obj,
    gr      = .clogit_grad,
    Y = Y, Z = Z, v = v, logP0 = logP0, logW0 = logW0,
    method  = "BFGS",
    control = list(maxit = con$maxit, reltol = con$reltol),
    hessian = FALSE   # analytic Hessian computed below
  )

  lambda <- fit$par
  beta   <- -lambda

  # ---- reconstruct primal solutions ----------------------------------------
  V  <- .clogit_utilities(lambda, Z)
  pe <- .clogit_compute_pe(V, v, logP0, logW0)
  P  <- pe$P; W <- pe$W; E <- pe$E

  # ---- SE and vcov (analytic Hessian — always PSD) -------------------------
  H_analytic <- .clogit_analytic_hessian(Z, P, W, v)
  vcov_mat <- tryCatch(
    solve(H_analytic),
    error = function(e) {
      warning("Analytic Hessian singular; SEs set to NA.")
      matrix(NA_real_, K, K)
    }
  )
  se_beta <- sqrt(pmax(diag(vcov_mat), 0))

  # ---- names ---------------------------------------------------------------
  names(lambda) <- names(beta) <- names(se_beta) <- d$names
  dimnames(vcov_mat) <- list(d$names, d$names)
  colnames(P) <- colnames(E) <- d$alt_levels

  # ---- entropy measures ----------------------------------------------------
  Plog  <- ifelse(P > 0, P * log(P), 0.0)
  S_p   <- -sum(Plog) / (log(J) * N)

  Wf    <- as.numeric(W)
  Wlog  <- ifelse(Wf > 0, Wf * log(Wf), 0.0)
  S_w   <- -sum(Wlog) / (log(M) * N * J)

  KL_p  <- sum(ifelse(P > 0, P * (log(P) - logP0), 0.0))
  KL_w  <- sum(ifelse(Wf > 0, Wf * (log(Wf) - as.numeric(logW0)), 0.0))

  misses <- sum(max.col(P, ties.method = "first") != max.col(Y))

  H_signal <- -rowSums(Plog)          # N-vector (CLAUDE.md required)
  H_error  <- -sum(Wlog)              # scalar total error entropy

  structure(
    list(
      call          = mc,
      formula       = formula,
      beta          = beta,
      lambda        = lambda,
      lambda_hat    = lambda,          # CLAUDE.md alias
      se_beta       = se_beta,
      se_lambda     = se_beta,         # identical
      vcov          = vcov_mat,
      p_hat         = P,
      w_hat         = W,
      fitted.values = P,
      e             = E,
      y_mat         = Y,
      Z             = Z,
      P0            = P0,
      W0            = W0,
      v             = v,
      M             = M,
      loglik        = -fit$value,      # conventional: larger = better
      objective     = fit$value,       # raw minimised value
      H_signal      = H_signal,
      H_error       = H_error,
      S_p           = S_p,
      S_w           = S_w,
      S_total       = S_p + S_w,
      S             = S_p + S_w,       # CLAUDE.md alias
      KL_p          = KL_p,
      KL_w          = KL_w,
      misses        = misses,
      hessian       = H_analytic,
      converged     = (fit$convergence == 0L),
      iterations    = fit$counts,
      method        = "dual",
      N = N, J = J, K = K, M = M,
      alt_levels    = d$alt_levels
    ),
    class = c("infometrics_clogit", "infometrics")
  )
}


# ===========================================================================
#  S3 methods
# ===========================================================================

#' @export
print.infometrics_clogit <- function(x, digits = 4L, ...) {
  uniform_p <- max(abs(x$P0 - 1 / x$J)) < 1e-10
  uniform_w <- max(abs(x$W0 - 1 / x$M)) < 1e-10
  estimator <- if (uniform_p && uniform_w) "GME" else "GCE"

  cat(sprintf("\nConditional Logit (%s)\n", estimator))
  cat(strrep("-", 45L), "\n")
  cat(sprintf("  N = %d   J = %d   K = %d   M = %d\n", x$N, x$J, x$K, x$M))
  cat(sprintf("  alternatives : %s\n",
              paste(x$alt_levels, collapse = ", ")))
  cat(sprintf("  prior p0     : %s\n",
              if (uniform_p) "uniform (1/J)" else "non-uniform"))
  cat(sprintf("  prior w0     : %s\n",
              if (uniform_w) "uniform (1/M)" else "non-uniform"))
  cat(sprintf("  support v    : [%.4f, %.4f]  M = %d\n",
              min(x$v), max(x$v), x$M))
  cat(sprintf("  log-lik      : %.4f\n", x$loglik))
  cat(sprintf("  S(p_hat)     : %.4f   KL(p||p0) = %.4f\n",
              x$S_p, x$KL_p))
  cat(sprintf("  S(w_hat)     : %.4f   KL(w||w0) = %.4f\n",
              x$S_w, x$KL_w))
  cat(sprintf("  misses       : %d / %d (%.1f%%)\n",
              x$misses, x$N, 100 * x$misses / x$N))
  cat(sprintf("  convergence  : %s\n", if (x$converged) "yes" else "no"))

  cat("\nCoefficients (beta = -lambda):\n")
  tab <- data.frame(
    Estimate = round(x$beta,     digits),
    Std.Err  = round(x$se_beta,  digits),
    t.stat   = round(x$beta / x$se_beta, digits),
    check.names = FALSE
  )
  print(tab)
  invisible(x)
}

#' @export
summary.infometrics_clogit <- function(object, digits = 4L, ...) {
  tab <- data.frame(
    Estimate = object$beta,
    Std.Err  = object$se_beta,
    t.stat   = object$beta / object$se_beta,
    check.names = FALSE
  )

  uniform_p <- max(abs(object$P0 - 1 / object$J)) < 1e-10
  uniform_w <- max(abs(object$W0 - 1 / object$M)) < 1e-10
  estimator <- if (uniform_p && uniform_w) "GME" else "GCE"

  cat(sprintf("\nConditional Logit (%s)\n", estimator))
  cat("Call: "); print(object$call); cat("\n")
  cat(sprintf("  N=%d  J=%d  K=%d  M=%d\n",
              object$N, object$J, object$K, object$M))
  cat(sprintf("  log-lik=%.4f  S(p)=%.4f  S(w)=%.4f  misses=%d/%d\n",
              object$loglik, object$S_p, object$S_w,
              object$misses, object$N))
  cat("\nCoefficients:\n")
  print(round(tab, digits))
  invisible(object)
}

#' @export
coef.infometrics_clogit <- function(object, ...) {
  object$beta
}

#' @export
fitted.infometrics_clogit <- function(object, ...) {
  object$fitted.values
}

#' @export
residuals.infometrics_clogit <- function(object, ...) {
  object$y_mat - object$p_hat - object$e
}


# ===========================================================================
#  gce_clogitWrap() — stargazer-compatible wrapper
# ===========================================================================

#' Wrap a \code{gce_clogit} fit for use with \code{stargazer}
#'
#' @description
#' Converts a \code{\link{gce_clogit}} result into an object that
#' \pkg{stargazer} can consume directly.  Pass the wrapped object wherever
#' you would normally pass a model to \code{stargazer()}.
#'
#' @details
#' \pkg{stargazer} dispatches on \code{class(model)[1]}, so the wrapper must
#' present itself as \code{"lm"}.  Rather than subclassing, the wrapper
#' constructs a \emph{genuine} \code{lm} object on synthetic data engineered
#' so that stargazer extracts the correct GCE statistics:
#'
#' \strong{Construction.}  Let \eqn{H = V^{-1}} (the analytic Hessian) and
#' \eqn{C = \mathrm{chol}(H)} (upper-triangular Cholesky factor).  Set
#' \deqn{X = \begin{pmatrix} C \\ 0_{(N-K)\times K} \end{pmatrix}, \quad
#'       y = \begin{pmatrix} C\hat\beta \\ \mathbf{1}_{N-K} \end{pmatrix}.}
#' Then \eqn{X'X = H}, \eqn{(X'X)^{-1} = V}, and \eqn{\hat\beta_{\mathrm{lm}}
#' = \hat\beta}.  The residual sum of squares equals \eqn{N-K}, giving
#' \eqn{s^2 = 1} and \eqn{\mathrm{SE} = \sqrt{\mathrm{diag}(V)}}.
#' Sample size is exactly \eqn{N}.
#'
#' \strong{Suppressing OLS-specific rows.}  The synthetic R\eqn{^2} and
#' residual standard error are artefacts of the construction and should be
#' hidden.  Use
#' \preformatted{
#' stargazer(wrap,
#'           omit.stat = c("rsq", "adj.rsq", "ser"),
#'           add.lines = list(
#'             c("Log-likelihood",   round(wrap$.gce_loglik,    2)),
#'             c("McFadden R2",      round(wrap$.gce_pseudo_r2, 4))
#'           ))
#' }
#' The fields \code{.gce_loglik}, \code{.gce_N}, \code{.gce_K}, and
#' \code{.gce_pseudo_r2} are stored on the returned object for this purpose.
#'
#' @param fit An object of class \code{"infometrics_clogit"} returned by
#'   \code{\link{gce_clogit}}.
#'
#' @return A genuine \code{"lm"} object whose \code{coefficients},
#'   standard errors, t-statistics, p-values, and \code{nobs} match the
#'   GCE estimates exactly.  Additional fields \code{.gce_loglik},
#'   \code{.gce_pseudo_r2}, \code{.gce_N}, and \code{.gce_K} hold the
#'   real model-fit statistics for use in \code{add.lines}.
#'
#' @seealso \code{\link{gce_clogit}}
#'
#' @examples
#' set.seed(42)
#' alts <- c("beach", "boat", "charter", "pier"); J <- length(alts); N <- 30L
#' df <- data.frame(
#'   id     = rep(seq_len(N), each = J),
#'   alt    = rep(alts, times = N),
#'   price  = round(runif(N * J, 50, 300)),
#'   catch  = round(runif(N * J, 0.1, 1.5), 2),
#'   income = rep(round(runif(N, 20, 80)), each = J)
#' )
#' eta  <- with(df, -0.01 * price + 0.5 * catch)
#' prob <- tapply(exp(eta), df$id, function(x) x / sum(x))
#' df$chosen <- unlist(lapply(seq_len(N), function(i)
#'   seq_len(J) == sample(J, 1L, prob = prob[[i]])))
#' fit  <- gce_clogit(chosen ~ price + catch | income, data = df,
#'                    alt = "alt", id = "id")
#' wrap <- gce_clogitWrap(fit)
#' coef(wrap)
#' sqrt(diag(stats::vcov(wrap)))   # == fit$se_beta
#'
#' @export
gce_clogitWrap <- function(fit) {
  if (!inherits(fit, "infometrics_clogit"))
    stop("fit must be an object of class \"infometrics_clogit\" ",
         "returned by gce_clogit().")

  K    <- fit$K
  N    <- fit$N
  beta <- fit$beta
  V    <- fit$vcov

  # ------------------------------------------------------------------
  # Build synthetic lm data (see @details for derivation):
  #   H = solve(V),  C = chol(H)
  #   X = rbind(C, 0_{(N-K) x K}),  y = c(C beta, 1_{N-K})
  # ------------------------------------------------------------------
  H <- tryCatch(
    solve(V),
    error = function(e) {
      warning("vcov is singular; falling back to diagonal precision.")
      diag(1 / pmax(diag(V), .Machine$double.eps))
    }
  )
  C <- tryCatch(
    chol(H),
    error = function(e) {
      warning("Precision matrix not positive definite; ",
              "using diagonal Cholesky approximation.")
      diag(sqrt(pmax(diag(H), 0)))
    }
  )

  n_extra <- max(N - K, 1L)         # at least 1 df so lm doesn't complain
  cnames  <- paste0(".gcez", seq_len(K))

  X_syn <- rbind(C, matrix(0.0, n_extra, K))
  y_syn <- c(drop(C %*% beta), rep(1.0, n_extra))
  colnames(X_syn) <- cnames

  syn_df      <- as.data.frame(X_syn)
  syn_df$.y   <- y_syn
  fml         <- stats::as.formula(
    paste(".y ~", paste(cnames, collapse = " + "), "- 1")
  )
  lm_fit <- stats::lm(fml, data = syn_df)

  # stargazer identifies lm objects by checking object$call[1] == "lm()".
  # stats::lm() stores the call as "stats::lm(...)" which fails that check.
  # Patch the call so stargazer recognises the object correctly.
  lm_fit$call[[1L]] <- as.name("lm")

  # Rename coefficients from .gcez1,.gcez2,... to GCE parameter names.
  # summary.lm uses names(coef(object)) for row labels, so this propagates
  # through to the stargazer table automatically.
  names(lm_fit$coefficients) <- names(beta)

  # ------------------------------------------------------------------
  # Attach real model-fit statistics as extra list elements.
  # Users can reference these in stargazer's add.lines argument.
  # ------------------------------------------------------------------
  loglik_cl          <- sum(log(pmax(fit$p_hat[fit$y_mat == 1L],
                                     .Machine$double.eps)))
  loglik_null        <- N * log(1 / fit$J)
  lm_fit$.gce_loglik    <- fit$loglik
  lm_fit$.gce_pseudo_r2 <- 1 - loglik_cl / loglik_null
  lm_fit$.gce_N         <- N
  lm_fit$.gce_K         <- K
  lm_fit$.gce_S_p       <- fit$S_p
  lm_fit$.gce_S_w       <- fit$S_w
  lm_fit$.gce_misses    <- fit$misses
  lm_fit$.gce_v         <- fit$v

  lm_fit
}
