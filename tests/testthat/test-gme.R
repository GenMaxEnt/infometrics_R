test_that("gme() converges and satisfies moment constraints on simulated data", {
  set.seed(1)
  n <- 50
  x <- rnorm(n)
  y <- 1 + 2 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)

  fit <- gme(y ~ x, data = df)

  expect_true(fit$converged)
  # moment constraints: y = y_hat + e_hat (CLAUDE.md tolerance: 1e-4 for dual solver)
  expect_lt(max(abs(y - fitted(fit) - residuals(fit))), 1e-4)
})

test_that("GME == GCE with uniform priors", {
  set.seed(2)
  n <- 40
  x <- rnorm(n)
  y <- 2 + 3 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)

  fit_gme <- gme(y ~ x, data = df)

  # explicit uniform priors — must give identical coefficients
  K  <- 2L; M  <- ncol(fit_gme$p_hat); J  <- ncol(fit_gme$w_hat)
  p0 <- matrix(1 / M, nrow = K, ncol = M)
  w0 <- matrix(1 / J, nrow = n, ncol = J)
  fit_gce <- gme(y ~ x, data = df,
                 Z = fit_gme$Z, p0 = p0,
                 V = fit_gme$V, w0 = w0)

  expect_lt(max(abs(coef(fit_gme) - coef(fit_gce))), 1e-8)
})

test_that("S3 methods return objects of correct type and length", {
  set.seed(3)
  n <- 40
  df <- data.frame(y = rnorm(n), x = rnorm(n))
  fit <- gme(y ~ x, data = df)

  expect_length(coef(fit), 2L)
  expect_length(fitted(fit), n)
  expect_length(residuals(fit), n)
  expect_equal(dim(vcov(fit)), c(2L, 2L))
  expect_output(print(fit))
  expect_output(summary(fit))
})

test_that("normalised entropy S_p is in [0, 1] for each coefficient", {
  set.seed(4)
  n <- 50
  df <- data.frame(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
  fit <- gme(y ~ x1 + x2, data = df)

  expect_true(all(fit$S_p >= 0 & fit$S_p <= 1))
  expect_true(fit$S_P >= 0 && fit$S_P <= 1)
  expect_true(fit$pseudo_R2 >= 0 && fit$pseudo_R2 <= 1)
})

test_that("p_hat rows and w_hat rows are valid probability vectors", {
  set.seed(5)
  n <- 40
  df <- data.frame(y = rnorm(n), x = rnorm(n))
  fit <- gme(y ~ x, data = df)

  expect_lt(max(abs(rowSums(fit$p_hat) - 1)), 1e-8)
  expect_lt(max(abs(rowSums(fit$w_hat) - 1)), 1e-8)
  expect_true(all(fit$p_hat > 0))
  expect_true(all(fit$w_hat > 0))
})

test_that("gme() input validation: nu out of range", {
  df <- data.frame(y = rnorm(20), x = rnorm(20))
  expect_error(gme(y ~ x, data = df, nu = 0),   "nu must be a single number in")
  expect_error(gme(y ~ x, data = df, nu = 1),   "nu must be a single number in")
  expect_error(gme(y ~ x, data = df, nu = -0.1),"nu must be a single number in")
  expect_error(gme(y ~ x, data = df, nu = 1.5), "nu must be a single number in")
})

test_that("gme() input validation: asymmetric V", {
  set.seed(6)
  n <- 30
  df <- data.frame(y = rnorm(n), x = rnorm(n))
  K  <- 2L
  Z  <- matrix(c(-2, -1, 0, 1, 2), nrow = K, ncol = 5, byrow = TRUE)
  V_bad <- c(-3, 0, 4)   # not symmetric: sum != 0
  w0 <- matrix(1/3, nrow = n, ncol = 3)
  expect_error(
    gme(y ~ x, data = df, Z = Z, V = V_bad, w0 = w0),
    "symmetric around zero"
  )
})

test_that("gme() input validation: p0 rows not summing to 1", {
  set.seed(7)
  n <- 30
  df <- data.frame(y = rnorm(n), x = rnorm(n))
  K  <- 2L; M <- 3L
  sp <- default_supports(df$y, model.matrix(~x, df))
  p0_bad <- matrix(0.5, nrow = K, ncol = M)   # rows sum to 1.5
  expect_error(
    gme(y ~ x, data = df, Z = sp$Z[, 1:M], p0 = p0_bad),
    "rows summing to 1"
  )
})

test_that("coefficients are named with model matrix column names", {
  set.seed(8)
  n <- 40
  df <- data.frame(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
  fit <- gme(y ~ x1 + x2, data = df)
  expect_equal(names(coef(fit)), c("(Intercept)", "x1", "x2"))
})
