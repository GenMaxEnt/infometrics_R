# Tests for inverse_noise() — noisy-moment GME/GCE via the formula interface,
# including de-duplicated output, signal+noise normalized entropy, Fano bounds,
# and standard errors (delta / sandwich / bootstrap).

make_noise_data <- function() {
  X <- matrix(c(1, 2, 3, 4, 5,
                5, 4, 3, 2, 1,
                2, 1, 5, 3, 4), nrow = 3, byrow = TRUE)   # T = 3, K = 5
  p_true <- c(.30, .25, .20, .15, .10)
  data.frame(y  = as.vector(X %*% p_true),
             s1 = X[, 1], s2 = X[, 2], s3 = X[, 3],
             s4 = X[, 4], s5 = X[, 5])
}

# Larger design (T = 10 > K = 3) for the SE tests (T = sample size).
make_reg_data <- function(seed = 123) {
  set.seed(seed)
  X <- matrix(runif(30), nrow = 10, ncol = 3)
  p_true <- c(0.5, 0.2, 0.3)
  y <- as.vector(X %*% p_true) + rnorm(nrow(X), 0, 0.3)
  data.frame(y = y, x1 = X[, 1], x2 = X[, 2], x3 = X[, 3])
}

f5 <- y ~ s1 + s2 + s3 + s4 + s5 - 1
f3 <- y ~ x1 + x2 + x3 - 1

test_that("converges and returns a well-formed object", {
  fit <- inverse_noise(f5, data = make_noise_data())

  expect_s3_class(fit, "inverse_noise")
  expect_s3_class(fit, "infometrics")
  expect_equal(fit$convergence, 0L)
  expect_equal(sum(coef(fit)), 1, tolerance = 1e-10)
  expect_equal(fit$method, "dual")
  expect_length(coef(fit), 5)
})

test_that("output is de-duplicated (canonical fields kept, aliases dropped)", {
  fit <- inverse_noise(f5, data = make_noise_data())
  nm  <- names(fit)
  # dropped duplicate/redundant fields
  expect_false(any(c("p", "lambda", "w", "e", "value", "converged") %in% nm))
  # kept canonical fields
  expect_true(all(c("p_hat", "lambda_hat", "w_hat", "residuals", "objective",
                    "convergence") %in% nm))
  # accessors read the kept fields
  expect_equal(unname(coef(fit)), unname(fit$p_hat))
  expect_equal(unname(residuals(fit)), unname(as.vector(fit$y - fit$fitted.values)))
})

test_that("noisy moment condition holds: residuals equal the recovered noise", {
  fit <- inverse_noise(f5, data = make_noise_data())
  # At the optimum grad = Xp + e - y = 0, so y - Xp = e = sum_j w_tj v_j.
  e_recon <- as.vector(fit$w_hat %*% fit$support)
  expect_lt(max(abs(residuals(fit) - e_recon)), 1e-6)
  expect_equal(unname(fitted(fit)), as.vector(fit$fitted.values),
               tolerance = 1e-12)
})

test_that("shrinking the noise support recovers the pure solution", {
  X <- rbind(c(1, 2, 3, 4, 5),
             c(1, 4, 9, 16, 25))
  p_true <- exp(as.vector(crossprod(X, c(0.2, -0.05))))
  p_true <- p_true / sum(p_true)
  dat <- data.frame(y  = as.vector(X %*% p_true),
                    s1 = X[, 1], s2 = X[, 2], s3 = X[, 3],
                    s4 = X[, 4], s5 = X[, 5])

  fit_pure  <- inverse_ce(f5, data = dat)
  fit_noise <- inverse_noise(f5, data = dat, v = c(-1e-5, 0, 1e-5))

  expect_equal(unname(coef(fit_noise)), unname(coef(fit_pure)),
               tolerance = 1e-3)
})

test_that("GME (uniform priors) differs from GCE (non-uniform signal prior)", {
  dat <- make_noise_data()
  fit_gme <- inverse_noise(f5, data = dat)                       # uniform p0
  fit_gce <- inverse_noise(f5, data = dat,
                           p0 = c(.50, .20, .15, .10, .05))      # skewed prior

  expect_false(isTRUE(all.equal(unname(coef(fit_gme)),
                                unname(coef(fit_gce)),
                                tolerance = 1e-3)))
  expect_equal(sum(coef(fit_gce)), 1, tolerance = 1e-10)
})

test_that("wider noise support shrinks the signal toward its prior", {
  dat <- make_noise_data()
  fit_narrow <- inverse_noise(f5, data = dat, v = c(-0.05, 0, 0.05))
  fit_wide   <- inverse_noise(f5, data = dat, v = c(-50, 0, 50))
  expect_gt(fit_wide$S, fit_narrow$S)
})

test_that("normalized entropy for signal AND noise are reported and in [0, 1]", {
  fit <- inverse_noise(f5, data = make_noise_data())
  expect_true(all(c("S", "S_p", "S_w", "H_signal", "H_error", "H_error0")
                  %in% names(fit)))
  expect_equal(fit$S, fit$S_p)
  expect_gte(fit$S, 0);   expect_lte(fit$S, 1)
  expect_gte(fit$S_w, 0); expect_lte(fit$S_w, 1)     # uniform w0 => S_w in [0,1]
  expect_equal(fit$S_w, fit$H_error / fit$H_error0, tolerance = 1e-10)
})

test_that("analytic Hessian is symmetric positive definite", {
  fit <- inverse_noise(f5, data = make_noise_data())
  H <- fit$hessian
  expect_equal(H, t(H), tolerance = 1e-10)
  expect_gt(min(eigen(H, symmetric = TRUE, only.values = TRUE)$values), 0)
})

# ---- standard errors -------------------------------------------------------

test_that("sandwich SEs are stored, finite, and positive by default", {
  fit <- inverse_noise(f3, data = make_reg_data())
  expect_equal(fit$se_method, "sandwich")
  expect_length(fit$se_p, 3)          # K
  expect_length(fit$se_lambda, 10)    # T
  expect_true(all(is.finite(fit$se_p) & fit$se_p > 0))
  expect_true(all(is.finite(fit$se_lambda) & fit$se_lambda > 0))
  expect_false(is.null(fit$vcov_lambda))
})

test_that("vcov(type=) returns the T x T covariance for each method", {
  fit <- inverse_noise(f3, data = make_reg_data())
  Vs <- vcov(fit, type = "sandwich")
  Vd <- vcov(fit, type = "delta")
  expect_equal(dim(Vs), c(10L, 10L))
  expect_equal(dim(Vd), c(10L, 10L))
  # delta is the raw inverse Hessian (the pre-change default)
  expect_equal(unname(Vd), unname(solve(fit$hessian)), tolerance = 1e-8)
  # sandwich is positive definite
  expect_gt(min(eigen(Vs, symmetric = TRUE, only.values = TRUE)$values), 0)
  set.seed(1)
  Vb <- vcov(fit, type = "bootstrap", B = 100)
  expect_equal(dim(Vb), c(10L, 10L))
})

test_that("the naive delta SEs overstate the sandwich SEs", {
  fit <- inverse_noise(f3, data = make_reg_data())
  sed <- summary(fit, se_method = "delta")$se_p
  ses <- summary(fit, se_method = "sandwich")$se_p
  expect_true(all(sed > ses))
  expect_gt(mean(sed / ses), 2)       # delta overstates substantially
})

test_that("residual bootstrap SEs are in the same ballpark as the sandwich", {
  fit <- inverse_noise(f3, data = make_reg_data())
  ses <- summary(fit, se_method = "sandwich")$se_p
  set.seed(2)
  seb <- summary(fit, se_method = "bootstrap", B = 200)$se_p
  expect_length(seb, 3)
  expect_true(all(seb > 0))
  r <- mean(seb / ses)
  expect_gt(r, 0.3); expect_lt(r, 3)  # loose band (both ~0.8-1.2x of MC truth)
})

test_that("se_method = 'none' skips SE computation", {
  fit <- inverse_noise(f3, data = make_reg_data(), se_method = "none")
  expect_equal(fit$se_method, "none")
  expect_null(fit$se_p)
  expect_null(fit$se_lambda)
  expect_null(fit$vcov_lambda)
})

# ---- Fano bounds -----------------------------------------------------------

test_that("fano_bounds returns a one-row frame with a valid Fano bound", {
  fit <- inverse_noise(f3, data = make_reg_data())
  fb  <- fano_bounds(fit)
  expect_s3_class(fb, "data.frame")
  expect_equal(nrow(fb), 1L)
  expect_true(all(c("p_max", "pe", "H", "S", "pe_lower") %in% names(fb)))
  K <- length(coef(fit))
  expect_equal(fb$pe, 1 - max(coef(fit)), tolerance = 1e-10)
  expect_gte(fb$pe, fb$pe_lower - 1e-9)
  # strong Fano: H(pe) + pe log(K-1) >= H(p)
  pe  <- fb$pe
  Hpe <- if (pe > 0 && pe < 1) -pe * log(pe) - (1 - pe) * log(1 - pe) else 0
  expect_gte(Hpe + pe * log(K - 1), fb$H - 1e-9)
})

# ---- input validation ------------------------------------------------------

test_that("input validation throws informative errors", {
  dat <- make_noise_data()

  expect_error(inverse_noise(f5, data = dat, nu = 0),   "strictly between 0 and 1")
  expect_error(inverse_noise(f5, data = dat, nu = 1),   "strictly between 0 and 1")
  expect_error(inverse_noise(f5, data = dat, nu = 1.5), "strictly between 0 and 1")
  expect_error(inverse_noise(f5, data = dat, p0 = rep(1 / 3, 3)),
               "must equal the number of states")
  expect_error(inverse_noise(f5, data = dat, p0 = c(.4, .3, .2, .1, 0)),
               "strictly positive")
  expect_error(inverse_noise(f5, data = dat, v = 1L),
               "must be >= 2")
  expect_error(
    inverse_noise(f5, data = dat, v = c(-1, 0, 1), w0 = matrix(1 / 2, 3, 2)),
    "must equal length\\(v\\)"
  )
  expect_error(
    inverse_noise(f5, data = dat, v = c(-1, 0, 1), w0 = matrix(1 / 3, 2, 3)),
    "must equal the number of moments"
  )
})

test_that("intercept in formula triggers a warning", {
  dat <- make_noise_data()
  expect_warning(
    inverse_noise(y ~ s1 + s2 + s3 + s4 + s5, data = dat),
    "intercept"
  )
})

test_that("the documented example runs without warnings", {
  dat <- make_reg_data(seed = 123)
  expect_no_warning(fit <- inverse_noise(f3, data = dat))
  expect_equal(fit$convergence, 0L)
})

test_that("scalar v is read as a support-point count", {
  fit <- inverse_noise(f5, data = make_noise_data(), v = 5)
  expect_length(fit$support, 5)
  expect_equal(mean(fit$support), 0, tolerance = 1e-10)   # symmetric
  expect_equal(dim(fit$w_hat), c(3L, 5L))
})

test_that("S3 methods run and return correctly shaped objects", {
  fit <- inverse_noise(f5, data = make_noise_data())

  expect_length(coef(fit), 5)
  expect_length(fitted(fit), 3)
  expect_length(residuals(fit), 3)
  expect_equal(dim(vcov(fit)), c(3L, 3L))                  # sandwich default
  expect_output(print(fit), "Estimated probabilities")
  s <- summary(fit)
  expect_s3_class(s, "summary.inverse_noise")
  expect_equal(s$se_method, "sandwich")
  expect_output(print(s), "estimated noise")
  expect_output(print(s), "Std. Error")                    # SE column
  expect_output(print(s), "t value")                       # t-stat column
  expect_output(print(s), "noise S\\(w\\)")                # noise entropy line
  expect_output(print(s), "Fano")                          # Fano line
})

test_that("summary honors an alternative se_method", {
  fit <- inverse_noise(f3, data = make_reg_data())
  s_d <- summary(fit, se_method = "delta")
  expect_equal(s_d$se_method, "delta")
  expect_output(print(s_d), "delta")
})

test_that("a healthy fit reports a near-zero first-order-condition residual", {
  fit <- inverse_noise(f3, data = make_reg_data())
  expect_true(is.numeric(fit$foc_residual))
  expect_lt(fit$foc_residual, 1e-4)
  # at the optimum the moment residual IS the estimated noise
  e_hat <- as.vector(fit$w_hat %*% fit$support)
  expect_lt(max(abs(as.vector(residuals(fit)) - e_hat)), 1e-4)
})

test_that("an infeasibly narrow noise support warns (unbounded dual)", {
  d <- make_reg_data()
  expect_warning(inverse_noise(f3, data = d, v = c(-1e-4, 0, 1e-4)),
                 "unbounded")
  fit <- suppressWarnings(inverse_noise(f3, data = d, v = c(-1e-4, 0, 1e-4)))
  expect_gt(fit$foc_residual, 1e-4)
})

test_that("default noise support widens with the number of moments", {
  set.seed(11)
  mk <- function(n) {
    X <- matrix(runif(3 * n), nrow = n, ncol = 3)
    data.frame(y = as.vector(X %*% c(0.5, 0.2, 0.3)) + rnorm(n, 0, 0.3),
               x1 = X[, 1], x2 = X[, 2], x3 = X[, 3])
  }
  f_small <- inverse_noise(f3, data = mk(10L),   se_method = "none")
  f_big   <- inverse_noise(f3, data = mk(2000L), se_method = "none")
  k_small <- max(f_small$support) / stats::sd(f_small$y)
  k_big   <- max(f_big$support)   / stats::sd(f_big$y)
  expect_equal(k_small, 3, tolerance = 1e-6)              # T=10 -> still 3
  expect_gt(k_big, k_small)
  expect_equal(k_big, sqrt(2 * log(2000)), tolerance = 1e-6)
})
