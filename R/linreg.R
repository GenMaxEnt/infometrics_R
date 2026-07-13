# linreg.R
# Generalized Cross-Entropy linear regression with an lm-style formula
# interface. Fits y = X beta + e by GCE: each coefficient is reparameterized
# on a bounded signal support and each error on a noise support, then the
# concentrated dual is solved in the T Lagrange multipliers.
#
# NOTE: linreg() fits the SAME estimator as gme() (GCE regression). It is a
# deliberate standalone sibling of inverse_ce()/inverse_noise() that adds
# lm-style conveniences (predict(), printCoefmat summary, r.squared) and keeps
# its own argument naming (lowercase `v` for the noise support). See CLAUDE.md
# "linreg vs gme" for the relationship.
#
# References:
#   Golan, A. (2008). Information and Entropy Econometrics -- A Review and
#     Synthesis. Foundations and Trends in Econometrics, 2(1-2), 1-145.
#     Section 6.1, p. 96 (asymptotic variance).

# ---- internal solver engine -------------------------------------------------
#
# Pure numeric core of the GCE-regression dual, factored out so that both the
# main fit (linreg) and the restricted refits used by the entropy-ratio test
# (summary.linreg) drive exactly the same solver path. All inputs are FULLY
# resolved (defaulting / recycling / validation happen in linreg() before this
# is called). Returns only what those two callers need; vcov / Hessian / R^2 /
# entropy post-processing stays in linreg() so restricted refits don't pay for
# it. No convergence warning here -- K refits would otherwise spam warnings.

#' @keywords internal
.linreg_engine <- function(X, y, Z, p0, v, w0, nu, con, init = NULL) {
  Tn    <- nrow(X)
  logp0 <- log(p0); logw0 <- log(w0)

  .signal <- function(lambda) {
    L <- logp0 + Z * (as.vector(lambda %*% X) / (1 - nu))
    rmax <- apply(L, 1L, max); EL <- exp(L - rmax); rs <- rowSums(EL)
    p <- EL / rs
    list(logOmega = rmax + log(rs), p = p, beta = as.vector(rowSums(Z * p)))
  }
  .noise <- function(lambda) {
    a <- logw0 + outer(as.vector(lambda), v) / nu
    amax <- apply(a, 1L, max); ea <- exp(a - amax); rs <- rowSums(ea)
    w <- ea / rs
    list(logPsi = amax + log(rs), w = w, e = as.vector(w %*% v))
  }
  ## minimised convex dual: -y'lambda + (1-nu) sum logOmega + nu sum logPsi
  obj <- function(lambda) {
    sg <- .signal(lambda); ns <- .noise(lambda)
    -sum(y * lambda) + (1 - nu) * sum(sg$logOmega) + nu * sum(ns$logPsi)
  }
  grad <- function(lambda) {
    sg <- .signal(lambda); ns <- .noise(lambda)
    as.vector(X %*% sg$beta + ns$e - y)
  }

  fit <- stats::optim(
    par     = if (is.null(init)) rep(0, Tn) else init,
    fn      = obj,
    gr      = grad,
    method  = "BFGS",
    control = list(maxit = con$maxit, reltol = con$reltol)
  )
  sg <- .signal(fit$par); ns <- .noise(fit$par)
  list(lambda = fit$par, p = sg$p, w = ns$w, beta = sg$beta, e = ns$e,
       value = fit$value, convergence = fit$convergence)
}

# Total joint entropy H*(fit) = sum H(p) + sum H(w) (natural log), over all
# entries of p (K x M) and w (T x J). Softmax outputs are normally strictly
# positive, but in a degenerate restricted refit (e.g. the overall test forces
# the intercept to 0 and the zero-centred noise support cannot span y, so the
# dual is unbounded and lambda diverges) the softmax underflows to exact zeros.
# Apply the standard 0*log0 = 0 convention by dropping zeros, exactly as
# shannon_entropy() and linreg()'s internal Hfun() do -- otherwise 0*log(0)
# yields NaN and the entropy-ratio statistic silently becomes NaN.

#' @keywords internal
.linreg_Hstar <- function(p, w) {
  pp <- p[p > 0]; ww <- w[w > 0]
  sum(-pp * log(pp)) + sum(-ww * log(ww))
}

# ---- estimator --------------------------------------------------------------

#' Generalized Cross-Entropy Linear Regression (Formula Interface)
#'
#' Fits the linear model \eqn{y = X\beta + e} by Generalized Cross-Entropy
#' (Golan 2008, Section 6.1, Eq. 6.5). Each coefficient is reparameterized on
#' a bounded support, \eqn{\beta_k = \sum_m z_{km} p_{km}}, and each error on a
#' noise support, \eqn{e_t = \sum_j v_j w_{tj}}; the signal and noise entropies
#' (relative to priors \code{p0}, \code{w0}) are maximized subject to the data
#' constraints, with weight \eqn{\nu \in (0,1)} on the noise. Estimation uses
#' the concentrated dual in the T Lagrange multipliers.
#'
#' As the supports widen the GCE estimates approach OLS; narrower supports
#' shrink \eqn{\beta} toward the support centers.
#'
#' @details
#' This is the same estimator family as \code{\link{gme}} (GCE regression),
#' exposed with an lm-style interface that adds a \code{\link{predict}} method,
#' an \code{lm}-style coefficient table, and \code{r.squared}. The dual is
#' minimised here (convex) to match the package convention used by
#' \code{\link{gme}}, \code{\link{inverse_ce}} and \code{\link{inverse_noise}};
#' this is equivalent to maximising the joint signal-plus-noise entropy. Uniform
#' priors (\code{p0 = NULL}, \code{w0 = NULL}) give GME; user priors give GCE.
#'
#' \code{summary()} reports a per-coefficient Entropy-Ratio (ER) test of
#' \eqn{H_0\!:\beta_k = 0} (Golan 2008, Sec. 6.4/6.6) in place of a Wald
#' \emph{z} test. To do so it refits the model \eqn{K+1} times (one per
#' coefficient plus one joint test of all coefficients = 0), each warm-started from the fitted
#' multipliers; the resolved \code{Z}, \code{p0}, \code{v}, \code{w0},
#' \code{X}, \code{y} are stored on the object to make these refits faithful.
#'
#' @param formula,data,subset,na.action Standard model-frame arguments, as in
#'   \code{\link[stats]{lm}}. The intercept is kept by default.
#' @param Z Coefficient support. Either a length-M numeric vector (a common
#'   support recycled across all coefficients) or a K-by-M matrix whose rows
#'   correspond, in order, to the columns of the model matrix (so row 1 is the
#'   intercept's support when an intercept is present). Default: a 5-point
#'   symmetric support spanning +/- a multiple of \code{sd(y)}.
#' @param nu Entropy weight in \eqn{(0,1)} on the noise term; default 0.5
#'   (equal-weight GME/GCE).
#' @param p0 Prior signal probabilities, a K-by-M matrix matching \code{Z};
#'   strictly positive rows summing to 1. Default uniform.
#' @param v,w0 Noise support and prior, as in \code{\link{inverse_noise}}:
#'   \code{v} is \code{NULL} (3-point default), a whole-number support count,
#'   or an explicit vector; \code{w0} is a strictly-positive T-by-J matrix with
#'   rows summing to 1 (default uniform).
#' @param control A named list merged over the defaults and passed to
#'   \code{\link[stats]{optim}} (BFGS): \code{maxit} (default 500) and
#'   \code{reltol} (default 1e-12).
#' @param ... Currently unused.
#'
#' @return An object of class \code{c("linreg", "infometrics")}: a list with
#'   \code{coefficients} (\eqn{\beta}, named by the model-matrix columns),
#'   \code{p}/\code{p_hat} (K-by-M signal probabilities), \code{lambda}/
#'   \code{lambda_hat}, \code{w}/\code{w_hat} (T-by-J noise probabilities),
#'   \code{e}, \code{nu}, \code{fitted.values}, \code{residuals}, \code{vcov}
#'   (asymptotic Cov(\eqn{\beta}), Golan 2008 p. 96), \code{hessian} (positive
#'   definite dual Hessian), \code{r.squared}, \code{H_signal} (K-vector of
#'   per-coefficient signal entropies), \code{S} (signal normalized entropy),
#'   \code{objective}/\code{value}, \code{converged}, \code{convergence},
#'   \code{method}, the resolved inputs \code{Z}, \code{p0}, \code{v},
#'   \code{w0}, \code{X}, \code{y}, \code{control} (retained so \code{summary()}
#'   can refit for the ER test), plus \code{terms}, \code{model},
#'   \code{xlevels}, \code{call}.
#'
#' @references Golan, A. (2008). \emph{Information and Entropy Econometrics -
#'   A Review and Synthesis}. Foundations and Trends in Econometrics,
#'   2(1-2), 1-145. Section 6.1.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' d$y <- 2 - 1.5 * d$x1 + 0.8 * d$x2 + rnorm(100, sd = 0.7)
#' fit <- linreg(y ~ x1 + x2, data = d, Z = seq(-20, 20, length.out = 5))
#' summary(fit)
#'
#' @seealso \code{\link{gme}} (the matrix/regression sibling of this estimator),
#'   \code{\link{inverse_ce}}, \code{\link{inverse_noise}}.
#' @importFrom stats model.frame model.matrix model.response optim sd
#' @importFrom stats printCoefmat delete.response .getXlevels predict
#' @export
linreg <- function(formula, data, Z = NULL, nu = 0.5,
                   p0 = NULL, v = NULL, w0 = NULL,
                   subset, na.action, control = list(), ...) {

  ## ---- model frame (lm-style) --------------------------------------------
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y  <- as.vector(stats::model.response(mf, "numeric"))
  X  <- stats::model.matrix(mt, mf)
  Tn <- nrow(X); K <- ncol(X)
  cn <- colnames(X)

  if (anyNA(y) || anyNA(X)) stop("Missing values present after na.action.")

  if (length(nu) != 1L || !is.finite(nu) || nu <= 0 || nu >= 1)
    stop("'nu' must be a single value strictly between 0 and 1. Got nu = ", nu,
         ".")

  ## ---- coefficient support Z and prior p0 (K x M) ------------------------
  sdY  <- stats::sd(y)
  span <- if (is.finite(sdY) && sdY > 0) 3 * sdY else 3
  if (is.null(Z)) Z <- seq(-span, span, length.out = 5L)
  if (is.null(dim(Z)))
    Z <- matrix(Z, nrow = K, ncol = length(Z), byrow = TRUE)   # shared support
  Z <- as.matrix(Z)
  if (nrow(Z) != K)
    stop("nrow(Z) (", nrow(Z), ") must equal the number of model-matrix ",
         "columns (", K, "). Rows of Z follow the columns: ",
         paste(cn, collapse = ", "), ".")
  M <- ncol(Z)
  rownames(Z) <- cn

  if (is.null(p0)) p0 <- matrix(1 / M, K, M)
  p0 <- as.matrix(p0)
  if (!all(dim(p0) == c(K, M)))
    stop("dim(p0) must match dim(Z) = ", K, " x ", M, ".")
  if (any(p0 <= 0)) stop("All signal prior masses 'p0' must be strictly positive.")
  rp <- rowSums(p0)
  if (any(abs(rp - 1) > 1e-8)) { warning("Rows of 'p0' renormalized."); p0 <- p0 / rp }

  ## ---- noise support v and prior w0 (T x J) ------------------------------
  is_count <- function(z) length(z) == 1L && is.finite(z) && abs(z - round(z)) < 1e-8
  if (is.null(v)) {
    J <- if (is.null(w0)) 3L else ncol(w0)
    v <- seq(-span, span, length.out = J)
  } else if (is_count(v)) {
    J <- as.integer(round(v))
    if (J < 2L) stop("Scalar 'v' is read as a support-point count and must be >= 2.")
    v <- seq(-span, span, length.out = J)
  } else { v <- as.vector(v); J <- length(v); if (J < 2L) stop("'v' needs >= 2 points.") }
  if (is.null(w0)) {
    w0 <- matrix(1 / J, Tn, J)
  } else {
    w0 <- as.matrix(w0)
    if (ncol(w0) != J) stop("ncol(w0) (", ncol(w0), ") must equal length(v) (", J, ").")
    if (nrow(w0) != Tn) stop("nrow(w0) (", nrow(w0), ") must equal n obs (", Tn, ").")
    if (any(w0 <= 0)) stop("All noise prior masses 'w0' must be strictly positive.")
    rw <- rowSums(w0)
    if (any(abs(rw - 1) > 1e-8)) { warning("Rows of 'w0' renormalized."); w0 <- w0 / rw }
  }
  if (abs(mean(v)) > 1e-8 * max(1, span))
    warning("Noise support 'v' is not centered at 0; noise will not be mean-zero.")

  ## ---- optimize via the shared engine ------------------------------------
  con <- list(maxit = 500L, reltol = 1e-12)
  con[names(control)] <- control
  eng <- .linreg_engine(X, y, Z, p0, v, w0, nu, con)
  if (eng$convergence != 0L)
    warning("optim did not converge (code ", eng$convergence, ").")

  lambda <- eng$lambda
  p <- eng$p; beta <- eng$beta; names(beta) <- cn; rownames(p) <- cn
  w <- eng$w; e <- eng$e
  yhat <- as.vector(X %*% beta); resid <- y - yhat

  ## ---- analytic dual Hessian (positive definite) ------------------------
  ## H = (1/(1-nu)) X diag(Var_p(Z)) X'  +  (1/nu) diag(Var_w(v))
  ## (Hessian of the convex minimised objective; PD thanks to the noise term.)
  varp <- rowSums(Z^2 * p) - beta^2          # length K: Var of Z under each p_k
  varw <- as.vector(w %*% (v^2)) - e^2       # length n: Var of v under each w_t
  H    <- (1 / (1 - nu)) * (X %*% (varp * t(X))) + (1 / nu) * diag(varw, Tn)

  ## ---- asymptotic Cov(beta): Golan (2008), p. 96 ------------------------
  ##   Var(beta) = [ sigma2(beta) / varpi2(beta) ] (X'X)^{-1}
  ##   sigma2(beta) = (1/n) sum_i lambda_i^2
  ##   varpi2(beta) = [ (1/n) sum_i ( sum_j v_j^2 w_ij
  ##                                  - (sum_j v_j w_ij)^2 ) ]^{-2}
  ## (book T = code n; common support v_j). The multiplier is taken on the
  ## standard-GME scale lambda_i/nu, since this dual divides the noise
  ## exponent by nu; at nu = 0.5 this is exactly the equal-weight GME for
  ## which the formula is derived, and Var(beta) -> OLS as supports widen.
  sigma2 <- mean((lambda / nu)^2)
  varpi2 <- mean(varw)^(-2)
  Vbeta  <- (sigma2 / varpi2) * solve(crossprod(X))
  dimnames(Vbeta) <- list(cn, cn)

  Hfun <- function(pr) { pr <- pr[pr > 0]; -sum(pr * log(pr)) }
  H_sig_k <- apply(p, 1L, Hfun)                            # K-vector
  S  <- mean(H_sig_k) / Hfun(rep(1 / M, M))               # avg signal norm. entropy
  sst <- sum((y - mean(y))^2); r2 <- 1 - sum(resid^2) / sst

  structure(
    list(
      coefficients  = beta,
      p             = p,
      p_hat         = p,                       # CLAUDE.md canonical alias
      lambda        = lambda,
      lambda_hat    = lambda,                  # CLAUDE.md canonical alias
      w             = w,
      w_hat         = w,                       # CLAUDE.md canonical alias
      e             = e,
      nu            = nu,
      fitted.values = yhat,
      residuals     = resid,
      vcov          = Vbeta,
      hessian       = H,
      r.squared     = r2,
      H_signal      = H_sig_k,                 # CLAUDE.md canonical field (K-vector)
      S             = S,
      objective     = eng$value,               # CLAUDE.md canonical field
      value         = eng$value,
      converged     = (eng$convergence == 0L), # CLAUDE.md canonical field
      convergence   = eng$convergence,
      method        = "dual",                  # CLAUDE.md canonical field
      ## resolved inputs retained so summary() can refit for the ER test
      Z             = Z,                        # K x M signal support
      p0            = p0,                       # K x M signal prior
      v             = v,                        # length-J noise support
      w0            = w0,                       # T x J noise prior
      X             = X,                        # T x K model matrix
      y             = y,                        # length-T response
      control       = con,                      # merged optim control
      terms         = mt,
      model         = mf,
      xlevels       = stats::.getXlevels(mt, mf),
      call          = cl
    ),
    class = c("linreg", "infometrics")
  )
}

## ----------------------------------------------------------------------------
## S3 methods (lm()-style)
## ----------------------------------------------------------------------------

#' @describeIn linreg Estimated coefficients \eqn{\beta}.
#' @param object,x A \code{linreg} object.
#' @export
coef.linreg <- function(object, ...) object$coefficients

#' @describeIn linreg Fitted values \eqn{X\beta}.
#' @export
fitted.linreg <- function(object, ...) object$fitted.values

#' @describeIn linreg Residuals \eqn{y - X\beta} (the estimated noise \code{e}).
#' @export
residuals.linreg <- function(object, ...) object$residuals

#' @describeIn linreg Asymptotic covariance of \eqn{\beta}, Golan (2008) p. 96
#'   (coincides with OLS as the supports widen).
#' @export
vcov.linreg <- function(object, ...) object$vcov

#' @describeIn linreg Predictions for \code{newdata} (or fitted values).
#' @param newdata Optional data frame of new observations.
#' @export
predict.linreg <- function(object, newdata, ...) {
  if (missing(newdata) || is.null(newdata)) return(object$fitted.values)
  tt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(tt, newdata, xlev = object$xlevels)
  Xn <- stats::model.matrix(tt, mf)
  as.vector(Xn %*% object$coefficients)
}

#' @describeIn linreg Compact display.
#' @param digits Significant digits to print.
#' @export
print.linreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  invisible(x)
}

#' @describeIn linreg lm()-style summary with an entropy-ratio coefficient
#'   table. For each coefficient \eqn{k} an Entropy-Ratio (ER) test of
#'   \eqn{H_0\!:\beta_k = 0} is computed by refitting with row \eqn{k} of the
#'   signal support collapsed to zero and forming
#'   \eqn{\mathrm{ER}_k = 2[H^*_{\mathrm{unrestricted}} -
#'   H^*_{\mathrm{restricted},k}]} with \eqn{H^* = \sum H(p) + \sum H(w)};
#'   under \eqn{H_0}, \eqn{\mathrm{ER}_k \sim \chi^2_1} (Golan 2008, Sec.
#'   6.4/6.6). An overall ER test of \eqn{H_0\!:} all coefficients (intercept
#'   and slopes) \eqn{= 0} is also reported. Note this performs \eqn{K+1}
#'   warm-started refits.
#' @importFrom stats pchisq
#' @export
summary.linreg <- function(object, ...) {
  beta <- object$coefficients
  cn   <- names(beta)
  K    <- length(beta)

  ## unrestricted joint entropy from the stored fit (no refit needed)
  Hstar_u <- .linreg_Hstar(object$p, object$w)

  ## per-coefficient ER test (all coefficients, including the intercept),
  ## warm-started from the unrestricted lambda for speed
  ER <- numeric(K); names(ER) <- cn
  for (k in seq_len(K)) {
    Zr <- object$Z; Zr[k, ] <- 0                # forces beta_k = sum(Zr[k,]*p_k) = 0
    eng_k <- .linreg_engine(object$X, object$y, Zr, object$p0, object$v,
                            object$w0, object$nu, object$control,
                            init = object$lambda)
    ER[k] <- 2 * (Hstar_u - .linreg_Hstar(eng_k$p, eng_k$w))
  }
  ER <- pmax(ER, 0)                             # clamp tiny negatives (-> p-value 1)
  pv <- stats::pchisq(ER, df = 1, lower.tail = FALSE)

  se  <- sqrt(diag(object$vcov))                # retained for reference + footer note
  tab <- cbind(Estimate = beta, `Std. Error` = se, ER = ER, `Pr(>Chi)` = pv)

  ## overall ER test: H0 all coefficients (intercept and slopes) = 0
  ## -> zero EVERY signal-support row, refit, df = K (Golan 2008, Sec. 6.4/6.6)
  Zr0 <- object$Z; Zr0[] <- 0
  eng0 <- .linreg_engine(object$X, object$y, Zr0, object$p0, object$v,
                         object$w0, object$nu, object$control,
                         init = object$lambda)
  ER_all <- max(0, 2 * (Hstar_u - .linreg_Hstar(eng0$p, eng0$w)))
  df_all <- K
  p_all  <- stats::pchisq(ER_all, df = df_all, lower.tail = FALSE)

  structure(list(call = object$call, coefficients = tab,
                 residuals = object$residuals, r.squared = object$r.squared,
                 nu = object$nu, S = object$S, convergence = object$convergence,
                 er_overall = c(ER = ER_all, df = df_all, p.value = p_all)),
            class = "summary.linreg")
}

#' @describeIn linreg Print method for the summary object.
#' @param signif.stars Logical; show significance stars (default from options).
#' @export
print.summary.linreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Residuals:\n"); print(summary(x$residuals, digits = digits))
  cat("\nCoefficients (entropy-ratio test, H0: beta_k = 0):\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.stars = signif.stars,
                      cs.ind = 1:2,        # Estimate, Std. Error
                      tst.ind = 3L,        # ER: test statistic
                      has.Pvalue = TRUE,   # last column (Pr(>Chi)) -> stars
                      na.print = "NA")
  cat("\nnu =", format(x$nu, digits = digits),
      "  R-squared:", format(x$r.squared, digits = digits),
      "  signal norm. entropy S:", format(x$S, digits = digits), "\n")
  eo <- x$er_overall
  if (!is.null(eo))
    cat(sprintf("Overall ER test (H0: all coefficients = 0): ER = %s, df = %d, Pr(>Chi) = %s\n",
                format(eo[["ER"]], digits = digits), as.integer(eo[["df"]]),
                format.pval(eo[["p.value"]], digits = digits)))
  cat("Convergence code:", x$convergence,
      if (x$convergence == 0L) "(converged)" else "(see ?optim)", "\n")
  cat("\nNote: ER is the entropy-ratio statistic for H0: beta_k = 0, computed as\n",
      "2*(Hstar_unrestricted - Hstar_restricted_k) where the restricted fit zeroes\n",
      "row k of the signal support and refits; under H0, ER ~ chi-squared(1) and\n",
      "Pr(>Chi) is the upper-tail p-value (clamped ER >= 0). Std. Error is the\n",
      "Golan (2008, p. 96) asymptotic SE, shown for reference.\n",
      sep = "")
  invisible(x)
}
