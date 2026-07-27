# Tests for linreg_iv() — stochastic-moments GME-IV (Golan 2008, pp. 89-91).

# Simulated endogenous-regressor data: xe correlated with the structural error u;
# z is a valid instrument (correlated with xe, uncorrelated with u).
sim_iv <- function(n = 400L, seed = 1, K_extra = 0L) {
  set.seed(seed)
  z  <- rnorm(n); u <- rnorm(n)
  xe <- 0.7 * z + u + rnorm(n)                 # endogenous
  W  <- if (K_extra > 0) matrix(rnorm(n * K_extra), n, K_extra) else NULL
  Xex <- cbind(1, xe, W)                       # design: intercept, endog, exog
  y  <- as.vector(Xex %*% c(1, 1.5, rep(0.5, K_extra))) + u
  IVm <- cbind(1, z, W)                        # exog regressors are own instruments
  list(y = y, X = Xex, IV = IVm, n = n, K = ncol(Xex))
}

wideZ <- function(K) matrix(c(-15, 0, 15), nrow = K, ncol = 3, byrow = TRUE)

test_that("converges and (nearly) satisfies the instrument moments", {
  d <- sim_iv()
  fit <- linreg_iv(d$y, d$X, d$IV, wideZ(d$K))

  expect_s3_class(fit, "linreg_iv")
  expect_s3_class(fit, "infometrics")
  expect_equal(fit$method, "dual")
  expect_true(fit$converged)
  expect_lt(fit$moment_resid, 1e-3)            # standardized moments ~ 0
})

test_that("recovers OLS when the instruments equal the regressors", {
  d <- sim_iv()
  fit <- linreg_iv(d$y, d$X, d$X, wideZ(d$K))   # IV = X
  b_ols <- solve(crossprod(d$X), crossprod(d$X, d$y))
  expect_lt(max(abs(unname(coef(fit)) - as.vector(b_ols))), 1e-2)
})

test_that("recovers just-identified IV (2SLS) and corrects OLS bias", {
  d <- sim_iv()
  fit <- linreg_iv(d$y, d$X, d$IV, wideZ(d$K))
  b_iv  <- solve(crossprod(d$IV, d$X), crossprod(d$IV, d$y))
  b_ols <- solve(crossprod(d$X),  crossprod(d$X,  d$y))
  expect_lt(max(abs(unname(coef(fit)) - as.vector(b_iv))), 1e-2)
  # the endogeneity correction: IV slope < biased OLS slope, near the truth 1.5
  expect_lt(abs(coef(fit)[2] - 1.5), abs(b_ols[2] - 1.5))
})

test_that("handles over-identified systems", {
  d <- sim_iv()
  IVover <- cbind(d$IV, d$IV[, 2]^2)           # extra instrument z^2
  fit <- linreg_iv(d$y, d$X, IVover, wideZ(d$K))
  expect_length(fit$lambda_hat, ncol(IVover))  # one multiplier per instrument
  expect_lt(fit$moment_resid, 1e-3)
  expect_true(fit$converged)
})

test_that("a narrow signal support shrinks beta toward the support center", {
  d <- sim_iv()
  f_wide   <- linreg_iv(d$y, d$X, d$IV, wideZ(d$K))
  f_narrow <- linreg_iv(d$y, d$X, d$IV, matrix(c(-1, 0, 1), d$K, 3, byrow = TRUE))
  expect_lt(abs(coef(f_narrow)[2]), abs(coef(f_wide)[2]))   # pulled toward 0
})

test_that("input validation throws informative errors", {
  d <- sim_iv(n = 80L)
  expect_error(linreg_iv(d$y, d$X, d$IV[, 1, drop = FALSE], wideZ(d$K)),
               "Under-identified")
  expect_error(linreg_iv(d$y, d$X[-1, ], d$IV, wideZ(d$K)), "must equal nrow")
  expect_error(linreg_iv(d$y, d$X, d$IV, wideZ(d$K), nu = 0), "in \\(0, 1\\)")
  expect_error(linreg_iv(d$y, d$X, d$IV, wideZ(d$K), nu = 1), "in \\(0, 1\\)")
  expect_error(linreg_iv(d$y, d$X, d$IV, matrix(0, d$K + 1, 3)),
               "nrow\\(Z\\)")
  expect_error(linreg_iv(d$y, d$X, d$IV, wideZ(d$K), p0 = matrix(1 / 2, d$K, 2)),
               "dim\\(p0\\)")
})

test_that("canonical aliases and S3 methods are present and consistent", {
  d <- sim_iv(n = 150L)
  fit <- linreg_iv(d$y, d$X, d$IV, wideZ(d$K))

  expect_identical(fit$coefficients, fit$b_hat)
  expect_identical(coef(fit), fit$b_hat)
  expect_identical(fit$objective, fit$value)
  expect_equal(fit$converged, fit$convergence == 0L)
  expect_length(fit$H_signal, d$K)
  expect_length(coef(fit), d$K)
  expect_length(fitted(fit), d$n)
  expect_length(residuals(fit), d$n)
  expect_equal(unname(fitted(fit) + residuals(fit)), d$y, tolerance = 1e-10)
  expect_output(print(fit), "GME-IV")
  expect_output(print(summary(fit)), "Coefficients")
})

test_that("standard errors: sandwich (default) finite, PD, K x K", {
  d <- sim_iv(n = 200L)
  fit <- linreg_iv(d$y, d$X, d$IV, wideZ(d$K))
  expect_identical(fit$se_method, "sandwich")
  expect_length(fit$se_beta, d$K)
  expect_true(all(is.finite(fit$se_beta)) && all(fit$se_beta > 0))
  V <- vcov(fit)
  expect_equal(dim(V), c(d$K, d$K))
  expect_equal(V, t(V), tolerance = 1e-10)                      # symmetric
  expect_gt(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values), 0)  # PD
  expect_equal(unname(fit$se_beta), unname(sqrt(diag(V))), tolerance = 1e-12)
})

test_that("delta and bootstrap SE methods run and are finite", {
  d <- sim_iv(n = 200L)
  fd <- linreg_iv(d$y, d$X, d$IV, wideZ(d$K), se_method = "delta")
  expect_identical(fd$se_method, "delta")
  expect_true(all(is.finite(fd$se_beta)) && all(fd$se_beta > 0))
  set.seed(7)
  fb <- linreg_iv(d$y, d$X, d$IV, wideZ(d$K), se_method = "bootstrap", boot = 40L)
  expect_identical(fb$se_method, "bootstrap")
  expect_true(all(is.finite(fb$se_beta)) && all(fb$se_beta > 0))
})

test_that("se_method = 'none' skips SEs; vcov() then recomputes", {
  d <- sim_iv(n = 200L)
  fn <- linreg_iv(d$y, d$X, d$IV, wideZ(d$K), se_method = "none")
  expect_null(fn$se_beta)
  expect_null(fn$vcov)
  expect_equal(dim(vcov(fn)), c(d$K, d$K))                      # recomputes (sandwich)
})

test_that("vcov(type=) switches the meat; delta differs from sandwich", {
  d <- sim_iv(n = 200L)
  fit <- linreg_iv(d$y, d$X, d$IV, wideZ(d$K))
  Vs <- vcov(fit, type = "sandwich"); Vd <- vcov(fit, type = "delta")
  expect_equal(dim(Vd), c(d$K, d$K))
  expect_false(isTRUE(all.equal(diag(Vs), diag(Vd))))
})

test_that("sandwich SE -> 2SLS robust SE as the support widens", {
  d <- sim_iv(n = 300L)
  fit <- suppressWarnings(
    linreg_iv(d$y, d$X, d$IV,
              matrix(c(-40, 0, 40), nrow = d$K, ncol = 3, byrow = TRUE)))
  b   <- solve(crossprod(d$IV, d$X), crossprod(d$IV, d$y))
  r   <- as.vector(d$y - d$X %*% b); br <- solve(crossprod(d$IV, d$X))
  se_2sls <- sqrt(diag(br %*% crossprod(d$IV, (r^2) * d$IV) %*% t(br)))
  expect_equal(unname(fit$se_beta), unname(se_2sls), tolerance = 0.03)
})

test_that("summary() prints an SE coefficient table with z and p-values", {
  d <- sim_iv(n = 150L)
  fit <- linreg_iv(d$y, d$X, d$IV, wideZ(d$K))
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Std. Error", out)))
  expect_true(any(grepl("z value", out)))
  expect_true(any(grepl("SE method: sandwich", out)))
})
