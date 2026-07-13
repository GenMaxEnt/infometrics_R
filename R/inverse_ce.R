# inverse_ce.R
# Pure (exact-moment) inverse problem via Maximum Entropy / Cross-Entropy,
# solved through the concentrated (dual) model. Formula interface sibling of
# me(): the response is the moment vector and each right-hand-side term is a
# state / support point.
#
# (inverse_pure() is retained as a deprecated alias -- see the bottom of this
# file. The cross-entropy method is what the name inverse_ce() now reflects:
# a uniform prior p0 gives Maximum Entropy, a non-uniform p0 gives Cross-Entropy.)
#
# References:
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145.
#     Section 4.2 (dual concentrated model, Eq. 4.4-4.7).

# ---- internal: rank-checked symmetric inverse -------------------------------

#' @keywords internal
.rank_checked_inverse <- function(H, tol = 1e-8) {
  ## Symmetric inverse of the information matrix, returning NULL when it is
  ## singular / rank-deficient (base R eigen; no MASS dependency). Used for
  ## Var(lambda) = I^{-1}: invert only if every eigenvalue exceeds tol * max,
  ## otherwise the null directions are unidentified and SEs must be NA.
  H  <- (H + t(H)) / 2                        # enforce exact symmetry
  ee <- eigen(H, symmetric = TRUE)
  vals <- ee$values
  if (length(vals) == 0L || max(vals) <= 0 || min(vals) <= tol * max(vals))
    return(NULL)
  ee$vectors %*% (t(ee$vectors) / vals)       # V diag(1/vals) V'
}

# ---- estimator --------------------------------------------------------------

#' Pure Inverse Problem via Cross-Entropy (Formula Interface)
#'
#' Estimates the K-dimensional distribution \code{p} closest (in
#' Kullback-Leibler divergence) to a prior \code{p0} subject to the T pure
#' moment constraints \eqn{y = X p}, by solving the unconstrained dual
#' (concentrated) model of Golan (2008, Section 4.2). This is the
#' \code{\link[stats]{lm}}-style formula sibling of \code{\link{me}}: both solve
#' the same maximum-entropy / cross-entropy problem and return identical
#' estimates for the same data. A uniform prior \code{p0} gives Maximum Entropy;
#' a non-uniform \code{p0} gives Cross-Entropy (hence the name).
#'
#' @details
#' The dual concentrated objective (Golan 2008, Eq. 4.4-4.5), minimised here to
#' match the package convention used by \code{\link{me}}, is
#' \deqn{\ell(\lambda) = -\sum_t \lambda_t y_t + \log \Omega(\lambda),}
#' with partition function
#' \eqn{\Omega(\lambda) = \sum_k p_{0,k} \exp(\sum_t \lambda_t x_{tk})}. Minimising
#' this convex objective is equivalent to maximising the entropy of \code{p}
#' subject to the moment constraints. With uniform \code{p0} this is Maximum
#' Entropy; with a non-uniform prior it is Cross-Entropy. The estimated
#' probabilities are recovered as
#' \deqn{\hat p_k = p_{0,k} \exp(\textstyle\sum_t \hat\lambda_t x_{tk}) /
#'   \Omega(\hat\lambda)} (Golan 2008, Eq. 4.3), and the analytic dual Hessian
#' is \eqn{\nabla^2 \ell = X \,\mathrm{diag}(p)\, X' - (Xp)(Xp)' =
#' \mathrm{Cov}_p(\text{moments})} (Eq. 4.7), which is positive semidefinite.
#'
#' @section Information-matrix standard errors:
#' The Hessian above \emph{is} the Fisher information matrix of the multipliers,
#' \eqn{I(\lambda) = \mathrm{Cov}_p(\text{moments})}, so its inverse gives their
#' covariance \eqn{\mathrm{Var}(\lambda) = I^{-1}(\lambda)} (Golan 2008, Eq.
#' ce-vcov-lambda; Cover & Thomas, 2006, Ch. 17). \code{inverse_ce()} reports the
#' resulting \code{se_lambda}, the delta-method probability SEs
#' \code{se_p} (\eqn{\mathrm{Var}(p) = J\,I^{-1}J'} with
#' \eqn{J_{kt} = \partial p_k/\partial\lambda_t = p_k(x_{tk}-E_p[x_t])}), and the
#' full \code{vcov_lambda} (via \code{\link{vcov}}). These are
#' \strong{information-matrix / curvature} quantities: they measure how well the
#' moments identify \eqn{\lambda} -- following the multiplier's reading as the
#' \emph{marginal information content of a moment} (a moment \eqn{t} with larger
#' \eqn{\mathrm{Var}_p(x_t)} carries more information, hence a smaller
#' \eqn{SE(\lambda_t)}; \eqn{\lambda_t\approx 0} means moment \eqn{t} adds little)
#' -- and are \emph{not} frequentist sampling standard errors, because a pure
#' inverse problem is deterministic. \eqn{I(\lambda)} is singular when a moment
#' row is constant/collinear or \eqn{T \ge K} (its rank is at most \eqn{K-1});
#' there the affected SEs are returned as \code{NA} rather than a pseudo-inverse's
#' misleadingly finite value. For sampling-based inference use the
#' stochastic-moment \code{\link{gme}} / \code{\link{inverse_noise}}.
#'
#' @section Formula orientation:
#' Unlike a regression, the response is the vector of \emph{moments} and each
#' right-hand-side term is a \emph{state}. The model matrix built from
#' \code{formula} and \code{data} is therefore the T-by-K matrix \eqn{X} with
#' rows = moment constraints (T) and columns = states / support points (K); the
#' data frame has one row per moment. The estimated coefficients (\code{coef})
#' are the K probabilities \code{p}, one per column. Because columns are states,
#' an intercept would add a spurious extra state; drop it with \code{- 1}. A
#' warning is issued if an intercept is present.
#'
#' @param formula A model formula. The response is the moment vector \code{y};
#'   each RHS term is a state. Usually written without an intercept, e.g.
#'   \code{y ~ x1 + x2 + x3 - 1}.
#' @param data A data frame (or environment) with one row per moment.
#' @param p0 Numeric vector of length K (number of model-matrix columns):
#'   strictly-positive prior probabilities for the Cross-Entropy formulation.
#'   Defaults to the uniform distribution, which makes CE reduce to ME. (Named
#'   \code{p0} to match the GME/GCE signal-prior convention; for \code{inverse_ce}
#'   it is a K-vector prior over the states, not a K-by-M support-point prior.)
#' @param subset Optional vector specifying a subset of moments (rows) to use,
#'   as in \code{\link[stats]{lm}}.
#' @param na.action A function indicating how to handle \code{NA}s in the model
#'   frame, as in \code{\link[stats]{lm}}.
#' @param control A named list of control parameters merged over the defaults
#'   and passed to \code{\link[stats]{optim}} (BFGS): \code{maxit} (default 500)
#'   and \code{reltol} (default 1e-10).
#' @param ... Currently unused.
#'
#' @return An object of class \code{c("inverse_ce", "infometrics")}: a list
#'   with components
#'   \describe{
#'     \item{\code{p_hat}}{Named K-vector of estimated probabilities.}
#'     \item{\code{lambda_hat}}{Named T-vector of Lagrange multipliers (dual
#'       variables).}
#'     \item{\code{fitted.values}}{Fitted moments \eqn{X p} (length T).}
#'     \item{\code{residuals}}{Moment residuals \eqn{y - X p} (length T).}
#'     \item{\code{hessian}}{Analytic T-by-T dual Hessian = information matrix
#'       \eqn{I(\lambda)=\mathrm{Cov}_p(\text{moments})} (Eq. 4.7), positive
#'       semidefinite.}
#'     \item{\code{vcov_lambda}}{\eqn{\mathrm{Var}(\lambda)=I^{-1}(\lambda)}
#'       (T-by-T), or \code{NULL} when \eqn{I} is singular.}
#'     \item{\code{se_lambda}}{Information-matrix SEs of \eqn{\lambda} (length T),
#'       or \code{NA} when \eqn{I} is singular.}
#'     \item{\code{se_p}}{Delta-method SEs of the probabilities \code{p}
#'       (length K), or \code{NA} when \eqn{I} is singular.}
#'     \item{\code{H_signal}}{Shannon entropy \eqn{H(\hat p)}.}
#'     \item{\code{S}}{Normalized entropy \eqn{H(\hat p)/H(p_0)} in [0, 1]
#'       (pseudo-R2 = 1 - S).}
#'     \item{\code{objective}}{Dual objective at the optimum.}
#'     \item{\code{p0}}{The prior distribution used.}
#'     \item{\code{convergence}}{Raw \code{optim} convergence code.}
#'     \item{\code{method}}{Solver used (\code{"dual"}).}
#'     \item{\code{terms}, \code{model}, \code{call}}{Standard model components.}
#'   }
#'
#' @references Golan, A. (2008). \emph{Information and Entropy Econometrics -
#'   A Review and Synthesis}. Foundations and Trends in Econometrics,
#'   2(1-2), 1-145. Section 4.2. Cover, T. M. & Thomas, J. A. (2006).
#'   \emph{Elements of Information Theory}, 2nd ed., Ch. 17.
#'
#' @seealso \code{\link{me}} for the matrix/vector interface to the same
#'   ME/CE estimator; \code{\link{gme}} and \code{\link{inverse_noise}} for the
#'   stochastic-moment (GME/GCE) estimators that report sampling standard errors;
#'   \code{\link{fano_bounds}} for the Fano error bound on the recovered \code{p}.
#'   \code{inverse_pure()} is a deprecated alias for \code{inverse_ce()}.
#'
#' @examples
#' # Recover a 5-state distribution from 2 exact moment constraints.
#' # The data frame has one row per moment (T = 2) and one column per
#' # state (K = 5). Because p_true is the maximum-entropy distribution for
#' # its own moments, it is recovered exactly.
#' X <- rbind(c(1, 2, 3, 4, 5),
#'            c(1, 4, 9, 16, 25))           # 2 moments x 5 states
#' p_true <- exp(as.vector(crossprod(X, c(0.2, -0.05))))
#' p_true <- p_true / sum(p_true)
#' dat <- data.frame(y  = as.vector(X %*% p_true),
#'                   s1 = X[, 1], s2 = X[, 2], s3 = X[, 3],
#'                   s4 = X[, 4], s5 = X[, 5])
#' fit <- inverse_ce(y ~ s1 + s2 + s3 + s4 + s5 - 1, data = dat)
#' fit
#' coef(fit)          # ~ p_true
#' summary(fit)       # info-matrix SEs of lambda / p + a Fano line
#' vcov(fit)          # Var(lambda) = I^{-1}
#' fano_bounds(fit)   # Fano error bound for the recovered p
#'
#' @importFrom stats model.frame model.matrix model.response optim vcov
#' @export
inverse_ce <- function(formula, data, p0 = NULL,
                       subset, na.action,
                       control = list(), ...) {

  ## ---- build the model frame the lm() way --------------------------------
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- stats::model.response(mf, "numeric")
  X <- stats::model.matrix(mt, mf)

  if (attr(mt, "intercept") == 1L)
    warning("Formula includes an intercept; it becomes an extra 'state' with ",
            "constant moment 1. Use '- 1' to drop it if unintended.")

  ## ---- coerce & validate -------------------------------------------------
  y  <- as.vector(y)
  Tn <- nrow(X)                 # number of moment constraints (Golan's T)
  K  <- ncol(X)                 # number of states / support points (K)

  if (!is.numeric(y))
    stop("The response (moments 'y') must be numeric.")
  if (anyNA(y) || anyNA(X))
    stop("Missing values present after na.action; clean 'data' or set ",
         "na.action.")
  if (length(y) != Tn)
    stop("length(y) (", length(y), ") must equal nrow(model matrix) (", Tn,
         ").")

  if (is.null(p0)) {
    p0 <- rep(1 / K, K)         # uniform prior  ==>  CE reduces to ME
  }
  p0 <- as.vector(p0)
  if (length(p0) != K)
    stop("length(p0) (", length(p0), ") must equal the number of states ",
         "(model-matrix columns = ", K, ").")
  if (any(p0 <= 0))
    stop("All prior masses 'p0' must be strictly positive (CE uses log(p/p0)).")
  if (abs(sum(p0) - 1) > 1e-8) {
    warning("Prior 'p0' did not sum to 1; renormalizing.")
    p0 <- p0 / sum(p0)
  }
  if (K <= Tn + 1L)
    warning("K = ", K, " is not greater than T + 1 = ", Tn + 1L,
            ": the problem is not under-determined and the dual Hessian is ",
            "rank-deficient (lambda may be non-unique).")

  logp0 <- log(p0)

  ## ---- control defaults --------------------------------------------------
  con <- list(maxit = 500L, reltol = 1e-10)
  con[names(control)] <- control

  ## ---- stable partition function & probabilities (log-sum-exp) -----------
  ## exponents z_k(lambda) = log p0_k + sum_t lambda_t x_tk   (length K)
  .z        <- function(lambda) logp0 + as.vector(lambda %*% X)
  .logOmega <- function(z) { mmax <- max(z); mmax + log(sum(exp(z - mmax))) }
  .probs    <- function(z) { e <- exp(z - max(z)); e / sum(e) }

  ## ---- dual objective (Eq. 4.4) and analytic gradient (Eq. 4.6) ----------
  ## Minimised convex objective: ell(lambda) = -lambda'y + log Omega(lambda).
  ## Equivalent to maximising entropy; matches me()'s sign convention.
  obj  <- function(lambda) { z <- .z(lambda); -sum(lambda * y) + .logOmega(z) }
  grad <- function(lambda) as.vector(X %*% .probs(.z(lambda)) - y)

  ## ---- optimise: dual is an unconstrained convex minimisation ------------
  fit <- stats::optim(
    par     = rep(0, Tn),
    fn      = obj,
    gr      = grad,
    method  = "BFGS",
    control = list(maxit = con$maxit, reltol = con$reltol)
  )
  if (fit$convergence != 0L)
    warning("optim did not converge (code ", fit$convergence, ").")

  lambda <- fit$par
  p      <- .probs(.z(lambda))
  names(p) <- colnames(X)
  names(lambda) <- if (is.null(rownames(X)))
    paste0("moment", seq_len(Tn)) else rownames(X)

  ## ---- analytic dual Hessian = information matrix (Eq. 4.7) ---------------
  ## I(lambda) = X diag(p) X' - (Xp)(Xp)' = Cov_p(moments)  (T x T, pos.
  ## semidefinite, the Hessian of the convex minimised dual objective).
  Ep      <- as.vector(X %*% p)                  # E_p[g], length T
  cov_mom <- X %*% (p * t(X)) - tcrossprod(Ep)   # = X diag(p) X' - Ep Ep'
  hessian <- cov_mom

  ## ---- information-matrix inference (Golan 2008, sec 4.2; Eq. 4.7) --------
  ## Var(lambda) = I(lambda)^{-1}. Curvature / identification quantities, NOT
  ## sampling SEs (the pure inverse problem is deterministic). I is singular
  ## when a moment row is constant/collinear or T >= K (rank <= K - 1); then the
  ## unidentified directions give NA rather than a pseudo-inverse's misleadingly
  ## finite value (solve-or-NA, not ginv).
  vcov_lambda <- .rank_checked_inverse(hessian)
  if (is.null(vcov_lambda)) {
    se_lambda <- rep(NA_real_, Tn)
    se_p      <- rep(NA_real_, K)
  } else {
    dimnames(vcov_lambda) <- list(names(lambda), names(lambda))
    se_lambda <- sqrt(pmax(diag(vcov_lambda), 0))
    ## delta method for the probabilities: Var(p) = J Var(lambda) J',
    ## J_{k,t} = dp_k/dlambda_t = p_k (x_tk - E_p[x_t]).
    J    <- p * (t(X) - matrix(Ep, K, Tn, byrow = TRUE))   # K x T
    varp <- J %*% vcov_lambda %*% t(J)                     # K x K
    se_p <- sqrt(pmax(diag(varp), 0))
  }
  names(se_lambda) <- names(lambda)
  names(se_p)      <- names(p)

  ## ---- normalized entropy  S = H(p) / H(p0)  (pseudo-R2 = 1 - S) ---------
  H_p  <- shannon_entropy(p)
  H_p0 <- shannon_entropy(p0)
  S    <- if (H_p0 > 0) H_p / H_p0 else NA_real_

  structure(
    list(
      p_hat         = p,                       # CLAUDE.md canonical field
      lambda_hat    = lambda,                  # CLAUDE.md canonical field
      fitted.values = Ep,
      residuals     = as.vector(y - Ep),
      hessian       = hessian,
      vcov_lambda   = vcov_lambda,
      se_lambda     = se_lambda,
      se_p          = se_p,
      H_signal      = H_p,                     # CLAUDE.md canonical field
      S             = S,
      objective     = fit$value,               # CLAUDE.md canonical field
      p0            = p0,                       # the ME/CE prior
      convergence   = fit$convergence,
      method        = "dual",                  # CLAUDE.md canonical field
      terms         = mt,
      model         = mf,
      call          = cl
    ),
    class = c("inverse_ce", "infometrics")
  )
}

## ----------------------------------------------------------------------------
## S3 methods (lm()-style)
## ----------------------------------------------------------------------------

#' @describeIn inverse_ce Estimated probabilities \code{p}.
#' @param object,x An \code{inverse_ce} object.
#' @export
coef.inverse_ce <- function(object, ...) object$p_hat

#' @describeIn inverse_ce Fitted moments \eqn{X p}.
#' @export
fitted.inverse_ce <- function(object, ...) object$fitted.values

#' @describeIn inverse_ce Moment-fitting residuals \eqn{y - X p}.
#' @export
residuals.inverse_ce <- function(object, ...) object$residuals

#' @describeIn inverse_ce Information-matrix covariance of the multipliers
#'   \eqn{\mathrm{Var}(\lambda)=I^{-1}(\lambda)} (T-by-T), or \code{NULL} when
#'   \eqn{I} is singular.
#' @export
vcov.inverse_ce <- function(object, ...) object$vcov_lambda

#' @describeIn inverse_ce Compact printed display.
#' @param digits Number of significant digits to print.
#' @export
print.inverse_ce <- function(x, digits = max(3L, getOption("digits") - 3L),
                             ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Estimated probabilities (p):\n")
  print.default(format(x$p_hat, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\nLagrange multipliers (lambda):\n")
  print.default(format(x$lambda_hat, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\nNormalized entropy S = ", format(x$S, digits = digits),
      "   (pseudo-R2 = ", format(1 - x$S, digits = digits), ")\n", sep = "")
  invisible(x)
}

#' @describeIn inverse_ce lm()-style summary (information-matrix SE tables + Fano).
#' @export
summary.inverse_ce <- function(object, ...) {
  K  <- length(object$p_hat)
  fb <- attr(.fano_row_bounds(matrix(object$p_hat, 1L, K), K), "overall")
  structure(
    list(call         = object$call,
         coefficients = object$p_hat,
         se_p         = object$se_p,
         lambda       = object$lambda_hat,
         se_lambda    = object$se_lambda,
         residuals    = object$residuals,
         S            = object$S,
         fano         = fb,
         singular     = is.null(object$vcov_lambda),
         convergence  = object$convergence),
    class = "summary.inverse_ce"
  )
}

#' @describeIn inverse_ce Print method for the summary object.
#' @export
print.summary.inverse_ce <- function(x,
                                      digits = max(3L, getOption("digits") - 3L),
                                      ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Moment-fitting residuals (y - Xp):\n")
  if (length(x$residuals) > 5L) {
    print(summary(x$residuals, digits = digits))
  } else {
    print.default(format(x$residuals, digits = digits), print.gap = 2L,
                  quote = FALSE)
  }

  cat("\nCoefficients (estimated probabilities p; delta-method SE from I^-1):\n")
  ptab <- cbind(Estimate = x$coefficients, `Std. Error` = x$se_p)
  print.default(round(ptab, digits), print.gap = 2L)

  cat("\nLagrange multipliers (lambda; info-matrix SE, Var(lambda) = I^-1):\n")
  z    <- x$lambda / x$se_lambda
  ltab <- cbind(Estimate = x$lambda, `Std. Error` = x$se_lambda,
                `z value` = z)
  print.default(round(ltab, digits), print.gap = 2L)

  cat("\nNormalized entropy S = ", format(x$S, digits = digits),
      "   pseudo-R2 = ", format(1 - x$S, digits = digits), "\n", sep = "")
  cat("Fano (sec 7.5): modal error pe = ", format(x$fano["mean_pe"], digits = digits),
      " >= bound ", format(x$fano["mean_pe_lower"], digits = digits), "\n", sep = "")
  cat("Convergence code:", x$convergence,
      if (x$convergence == 0L) "(converged)" else "(see ?optim)", "\n")

  cat("\nNote: the SEs above are information-matrix (curvature) quantities,\n",
      "Var(lambda) = I^-1 (Golan 2008 sec 4.2), measuring how well the moments\n",
      "identify lambda -- not frequentist sampling SEs (a pure inverse problem\n",
      "is deterministic). ",
      if (isTRUE(x$singular))
        "I(lambda) is singular here (constant/collinear moments\nor T >= K), so the SEs are NA. " else "",
      "For sampling-based inference use the\nstochastic-moment gme() / inverse_noise().\n",
      sep = "")
  invisible(x)
}

## ----------------------------------------------------------------------------
## Fano error bounds for the recovered distribution p  (Golan 2008, sec 7.5)
## ----------------------------------------------------------------------------

#' Fano Error Bounds for a Pure Inverse-Problem (Cross-Entropy) Fit
#'
#' The recovered \code{p} is a single distribution over the \eqn{K} states, so
#' Golan's (2008, sec 7.5) Fano bound applies to it directly and is the natural
#' companion to the normalized entropy \code{S}: reading \code{S} as an
#' irreducible prediction error, a modal classifier that guesses
#' \eqn{\arg\max_k p_k} errs with probability \eqn{pe = 1 - \max_k p_k}, and
#' Fano's inequality lower-bounds this by
#' \deqn{pe \ge S(p) - \log 2/\log K, \qquad S(p) = H(p)/\log K.}
#' The \code{S} used \emph{here} is the uniform-reference \eqn{H(p)/\log K} the
#' Fano bound requires -- distinct from \code{inverse_ce}'s prior-relative
#' reported \eqn{S = H(p)/H(p_0)}. This is an information-theoretic error bound on
#' prediction accuracy, \emph{not} a sampling standard error (for the
#' identification SEs of the multipliers see \code{\link{vcov.inverse_ce}} /
#' \code{summary}).
#'
#' @param object An \code{inverse_ce} object.
#' @param ... Unused.
#' @return A one-row data frame with columns \code{p_max}, \code{pe} (modal
#'   error), \code{H} (entropy, nats), \code{S} (uniform-normalized entropy), and
#'   \code{pe_lower} (Fano weak lower bound). An \code{"overall"} attribute holds
#'   \code{mean_pe}, \code{mean_pe_lower}, and \code{S_system} (here identical to
#'   the single row).
#' @references Golan, A. (2008). Information and Entropy Econometrics.
#'   \emph{Foundations and Trends in Econometrics}, \strong{2}(1-2), 1-145
#'   (sec 7.5); Fano, R. (1961). \emph{Transmission of Information}.
#' @seealso \code{\link{fano_bounds}}, \code{\link{inverse_ce}}
#' @export
fano_bounds.inverse_ce <- function(object, ...) {
  p <- object$p_hat
  .fano_row_bounds(matrix(p, 1L, length(p)), length(p))
}

## ----------------------------------------------------------------------------
## Deprecated alias
## ----------------------------------------------------------------------------

#' @describeIn inverse_ce Deprecated alias for \code{inverse_ce()}; forwards all
#'   arguments and returns an \code{inverse_ce} object (with a deprecation
#'   warning). The historical prior argument \code{q} is still accepted and
#'   mapped to \code{p0}.
#' @export
inverse_pure <- function(formula, data, p0 = NULL,
                         subset, na.action,
                         control = list(), ..., q = NULL) {
  .Deprecated("inverse_ce")
  cl <- match.call()
  cl[[1L]] <- quote(inverse_ce)
  ## map the historical prior name 'q' -> 'p0'
  if ("q" %in% names(cl)) {
    if ("p0" %in% names(cl))
      stop("Supply only one of 'p0' (preferred) or the deprecated 'q'.")
    names(cl)[names(cl) == "q"] <- "p0"
  }
  eval(cl, parent.frame())
}
