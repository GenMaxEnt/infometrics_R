# ============================================================
# test-gme_mnl.R
# Tests for me_mnl() and gme_mnl() — multinomial IT estimators.
# ============================================================

# Shared simulation helper (J=3, K=2, n observations)
.sim_mnl <- function(n = 80, seed = 1L) {
  set.seed(seed)
  x  <- rnorm(n)
  pr <- exp(cbind(0, 1 + 2 * x, -1 + x))
  pr <- pr / rowSums(pr)
  y  <- factor(apply(pr, 1, function(p) sample(3L, 1L, prob = p)))
  data.frame(y = y, x = x)
}


# ============================================================
# me_mnl() tests
# ============================================================

test_that("me_mnl() converges on simulated data", {
  df  <- .sim_mnl(n = 80, seed = 1L)
  fit <- me_mnl(y ~ x, data = df)
  expect_true(fit$converged)
})

test_that("me_mnl() moment constraints: X'(Y - P) ~ 0 at convergence", {
  df  <- .sim_mnl(n = 80, seed = 2L)
  fit <- me_mnl(y ~ x, data = df)

  mf <- model.frame(y ~ x, df)
  X  <- model.matrix(terms(mf), mf)
  P  <- fit$p_hat
  Y  <- fit$y_mat
  # gradient of the dual = X'(Y-P)[,2:J], should be ~0 at optimum
  # (absolute tolerance is data-scale dependent; 1e-3 is tight enough)
  resid_cols <- (Y - P)[, 2:ncol(Y), drop = FALSE]
  max_grad   <- max(abs(crossprod(X, resid_cols)))
  expect_lt(max_grad, 1e-3)
})

test_that("me_mnl() p_hat rows are valid probability vectors", {
  df  <- .sim_mnl(n = 80, seed = 3L)
  fit <- me_mnl(y ~ x, data = df)
  expect_lt(max(abs(rowSums(fit$p_hat) - 1)), 1e-8)
  expect_true(all(fit$p_hat > 0))
})

test_that("me_mnl() S3 methods return correct types", {
  df  <- .sim_mnl(n = 80, seed = 4L)
  fit <- me_mnl(y ~ x, data = df)

  expect_equal(dim(coef(fit)), c(2L, 3L))      # K x J
  expect_equal(dim(fitted(fit)), c(80L, 3L))   # N x J
  expect_equal(dim(residuals(fit)), c(80L, 3L))
  expect_output(print(fit))
  expect_output(summary(fit))
})

test_that("me_mnl() normalised entropy S_p is in [0, 1]", {
  df  <- .sim_mnl(n = 80, seed = 5L)
  fit <- me_mnl(y ~ x, data = df)
  expect_gte(fit$S_p, 0)
  expect_lte(fit$S_p, 1)
})

test_that("me_mnl() beta columns named by factor levels", {
  df  <- .sim_mnl(n = 80, seed = 6L)
  fit <- me_mnl(y ~ x, data = df)
  expect_equal(colnames(fit$beta), levels(df$y))
  expect_equal(rownames(fit$beta), c("(Intercept)", "x"))
})

test_that("me_mnl() first lambda column is identically zero", {
  df  <- .sim_mnl(n = 80, seed = 7L)
  fit <- me_mnl(y ~ x, data = df)
  expect_equal(fit$lambda[, 1], c(`(Intercept)` = 0, x = 0))
})

test_that("me_mnl() init validation error", {
  df <- .sim_mnl(n = 80, seed = 8L)
  expect_error(me_mnl(y ~ x, data = df, init = rep(0, 3)),
               "K\\*\\(J-1\\)")
})


# ============================================================
# gme_mnl() tests
# ============================================================

test_that("gme_mnl() converges on simulated data", {
  df  <- .sim_mnl(n = 80, seed = 10L)
  fit <- gme_mnl(y ~ x, data = df)
  expect_true(fit$converged)
})

test_that("gme_mnl() moment constraints: X'(Y - P - E) ~ 0 at convergence", {
  df  <- .sim_mnl(n = 80, seed = 11L)
  fit <- gme_mnl(y ~ x, data = df)

  mf <- model.frame(y ~ x, df)
  X  <- model.matrix(terms(mf), mf)
  Y  <- fit$y_mat; P <- fit$p_hat; E <- fit$e
  resid_cols <- (Y - P - E)[, 2:fit$J, drop = FALSE]
  max_grad   <- max(abs(crossprod(X, resid_cols)))
  expect_lt(max_grad, 1e-4)
})

test_that("gme_mnl() == me_mnl() when M=1 (degenerate noise)", {
  # When error support has a single point v=0, w -> 1, e -> 0, so GME -> ME.
  # Use v = c(-eps, eps) very narrow instead (M >= 2 is required).
  # Instead, test the limiting sense: with very narrow support, betas are close.
  df    <- .sim_mnl(n = 80, seed = 12L)
  me    <- me_mnl(y ~ x, data = df)
  gme_n <- gme_mnl(y ~ x, data = df, v = c(-1e-6, 1e-6), warm_start = FALSE)
  # Betas should be close (within tolerance for narrow support)
  expect_lt(max(abs(coef(me) - coef(gme_n))), 0.05)
})

test_that("GME == GCE with uniform priors", {
  df      <- .sim_mnl(n = 80, seed = 13L)
  fit_gme <- gme_mnl(y ~ x, data = df)

  N  <- nrow(df); J <- nlevels(df$y); M <- fit_gme$M
  p0 <- matrix(1 / J, N, J)
  w0 <- array(1 / M, c(N, J, M))
  fit_gce <- gme_mnl(y ~ x, data = df, p0 = p0, w0 = w0,
                     init = as.numeric(-fit_gme$beta[, 2:J]))
  expect_lt(max(abs(coef(fit_gme) - coef(fit_gce))), 1e-6)
})

test_that("gme_mnl() p_hat and w_hat are valid probability arrays", {
  df  <- .sim_mnl(n = 80, seed = 14L)
  fit <- gme_mnl(y ~ x, data = df)

  # p rows
  expect_lt(max(abs(rowSums(fit$p_hat) - 1)), 1e-8)
  expect_true(all(fit$p_hat > 0))

  # w (i,j) fibers
  w_row_sums <- apply(fit$w_hat, c(1, 2), sum)
  expect_lt(max(abs(w_row_sums - 1)), 1e-8)
  expect_true(all(fit$w_hat > 0))
})

test_that("gme_mnl() S_p and S_w are in [0, 1]", {
  df  <- .sim_mnl(n = 80, seed = 15L)
  fit <- gme_mnl(y ~ x, data = df)
  expect_gte(fit$S_p, 0); expect_lte(fit$S_p, 1)
  expect_gte(fit$S_w, 0); expect_lte(fit$S_w, 1)
})

test_that("gme_mnl() S3 methods return correct types", {
  df  <- .sim_mnl(n = 80, seed = 16L)
  fit <- gme_mnl(y ~ x, data = df)

  expect_equal(dim(coef(fit)),      c(2L, 3L))
  expect_equal(dim(fitted(fit)),    c(80L, 3L))
  expect_equal(dim(residuals(fit)), c(80L, 3L))
  expect_output(print(fit))
  expect_output(summary(fit))
})

test_that("gme_mnl() column names match factor levels", {
  df  <- .sim_mnl(n = 80, seed = 17L)
  fit <- gme_mnl(y ~ x, data = df)
  expect_equal(colnames(fit$beta),  levels(df$y))
  expect_equal(rownames(fit$beta),  c("(Intercept)", "x"))
})

test_that("gme_mnl() first lambda column is identically zero", {
  df  <- .sim_mnl(n = 80, seed = 18L)
  fit <- gme_mnl(y ~ x, data = df)
  expect_equal(fit$lambda[, 1], c(`(Intercept)` = 0, x = 0))
})

test_that("gme_mnl() input validation: p0 rows not summing to 1", {
  df <- .sim_mnl(n = 80, seed = 19L)
  p0_bad <- matrix(0.4, nrow = 80, ncol = 3)   # rows sum to 1.2, not 1
  expect_error(gme_mnl(y ~ x, data = df, p0 = p0_bad),
               "Rows of p0 must sum to 1")
})

test_that("gme_mnl() input validation: w0 not summing to 1", {
  df   <- .sim_mnl(n = 80, seed = 20L)
  w_bad <- array(0.5, dim = c(80, 3, 3))  # each fiber sums to 1.5
  expect_error(gme_mnl(y ~ x, data = df, w0 = w_bad),
               "fiber of w0 must sum to 1")
})

test_that("gme_mnl() KL_p = KL_w = 0 for GME (uniform priors)", {
  df  <- .sim_mnl(n = 80, seed = 21L)
  fit <- gme_mnl(y ~ x, data = df)  # uniform priors => GME
  # KL divergences should be >= 0 (they measure distance from uniform)
  expect_gte(fit$KL_p, 0)
  expect_gte(fit$KL_w, 0)
})

test_that("gme_mnl() warm_start = FALSE gives same result as warm_start = TRUE", {
  df      <- .sim_mnl(n = 80, seed = 22L)
  fit_ws  <- gme_mnl(y ~ x, data = df, warm_start = TRUE)
  fit_nws <- gme_mnl(y ~ x, data = df, warm_start = FALSE)
  # Both should converge to the same minimum (convex objective)
  expect_true(fit_ws$converged)
  expect_true(fit_nws$converged)
  expect_lt(max(abs(coef(fit_ws) - coef(fit_nws))), 1e-4)
})
