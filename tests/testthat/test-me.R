test_that("me: uniform solution when only normalization constraint", {
  # With no moment constraints (T=0 is not valid), but with one moment equal
  # to the prior mean, ME should recover near-uniform probabilities.
  K <- 6
  X <- matrix(seq_len(K), nrow = 1)   # support = 1:6, T = 1
  y <- sum((1:6) * rep(1/6, 6))       # mean of uniform = 3.5

  fit <- me(y = y, X = X)

  # Should recover (near) uniform distribution
  expect_equal(fit$p_hat, rep(1/6, K), tolerance = 1e-6)
  expect_equal(fit$S, 1, tolerance = 1e-6)
  expect_true(fit$converged)
})

test_that("me: higher observed mean shifts mass to higher outcomes", {
  X <- matrix(1:6, nrow = 1)

  fit_low  <- me(y = 2.5, X = X)
  fit_high <- me(y = 5.0, X = X)

  # Expected value should match constraint
  ev_low  <- sum(1:6 * fit_low$p_hat)
  ev_high <- sum(1:6 * fit_high$p_hat)
  expect_equal(ev_low,  2.5, tolerance = 1e-6)
  expect_equal(ev_high, 5.0, tolerance = 1e-6)

  # Higher mean -> more mass on higher outcomes
  expect_gt(fit_high$p_hat[6], fit_low$p_hat[6])
})

test_that("me: moment constraints satisfied at solution", {
  K <- 8
  X <- rbind(seq_len(K), seq_len(K)^2)   # two moments: mean and E[X^2]
  y <- c(4, 20)                            # target mean and E[X^2]

  fit <- me(y = y, X = X)

  fitted_moments <- as.numeric(X %*% fit$p_hat)
  expect_equal(fitted_moments, y, tolerance = 1e-4)
})

test_that("me: CE with uniform prior equals ME", {
  K <- 6
  X <- matrix(1:6, nrow = 1)
  y <- 4.5
  q <- rep(1/K, K)

  fit_me <- me(y = y, X = X)
  fit_ce <- me(y = y, X = X, q = q)

  expect_equal(fit_me$p_hat, fit_ce$p_hat, tolerance = 1e-8)
})

test_that("me: CE shifts solution toward non-uniform prior", {
  K <- 6
  X <- matrix(1:6, nrow = 1)
  y <- 3.5  # same as uniform mean

  # Prior that puts most weight on outcome 6
  q <- c(0.05, 0.05, 0.10, 0.10, 0.20, 0.50)

  fit_me <- me(y = y, X = X)            # ME: starts from uniform
  fit_ce <- me(y = y, X = X, q = q)    # CE: starts from skewed prior

  # CE should put more weight on outcome 6 than ME
  expect_gt(fit_ce$p_hat[6], fit_me$p_hat[6])
})

test_that("me: normalized entropy S is in [0, 1]", {
  X <- matrix(1:6, nrow = 1)
  for (y_val in c(2, 3.5, 5)) {
    fit <- me(y = y_val, X = X)
    expect_gte(fit$S, 0 - 1e-10)
    expect_lte(fit$S, 1 + 1e-10)
  }
})

test_that("me: S3 methods work without error", {
  X   <- matrix(1:6, nrow = 1)
  fit <- me(y = 4.0, X = X)

  expect_output(print(fit))
  expect_output(summary(fit))
  expect_length(coef(fit), 6)
  expect_length(fitted(fit), 1)
  expect_length(residuals(fit), 1)
})

test_that("me: residuals near zero at convergence", {
  K <- 8
  X <- rbind(seq_len(K), seq_len(K)^2)
  y <- c(4, 20)
  fit <- me(y = y, X = X)
  expect_equal(residuals(fit), c(0, 0), tolerance = 1e-6)
})

test_that("me: input validation — K must exceed T", {
  X <- matrix(1:3, nrow = 1)   # K = 3, T = 1, OK
  expect_no_error(me(y = 2, X = X))

  X_bad <- matrix(1:2, nrow = 2)   # K = 1, T = 2: K <= T, invalid
  expect_error(me(y = c(1, 2), X = X_bad), "under-determined")
})

test_that("me: input validation — dimension mismatch", {
  X <- matrix(1:6, nrow = 1)
  expect_error(me(y = c(1, 2), X = X), "nrow\\(X\\)")
})
