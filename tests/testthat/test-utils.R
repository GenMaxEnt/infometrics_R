test_that("normalize_data: max scaling gives values in [0, 1]", {
  x <- c(10, 30, 50, 20, 40)
  nx <- normalize_data(x)
  expect_lte(max(nx), 1 + 1e-12)
  expect_gte(min(nx), 0 - 1e-12)
  expect_equal(max(nx), 1, tolerance = 1e-12)
})

test_that("normalize_data: range scaling spans [0, 1]", {
  x <- c(5, 15, 25, 10)
  nx <- normalize_data(x, by = "range")
  expect_equal(min(nx), 0, tolerance = 1e-12)
  expect_equal(max(nx), 1, tolerance = 1e-12)
})

test_that("normalize_data: constant vector handled gracefully", {
  x <- c(3, 3, 3)
  expect_no_error(normalize_data(x))
  expect_no_error(normalize_data(x, by = "range"))
})

test_that("make_support: returns M equally spaced points", {
  s <- make_support(half_range = 3, M = 5)
  expect_length(s, 5)
  expect_equal(s[1], -3, tolerance = 1e-12)
  expect_equal(s[5],  3, tolerance = 1e-12)
  diffs <- diff(s)
  expect_equal(diffs, rep(diffs[1], length(diffs)), tolerance = 1e-12)
})

test_that("make_support: symmetric around center", {
  s <- make_support(half_range = 2, M = 5, center = 0)
  expect_equal(sum(s), 0, tolerance = 1e-12)
})

test_that("make_support: non-zero center", {
  s <- make_support(half_range = 2, M = 5, center = 10)
  expect_equal(s[1], 8,  tolerance = 1e-12)
  expect_equal(s[5], 12, tolerance = 1e-12)
})

test_that("make_support: M = 2 gives just endpoints", {
  s <- make_support(half_range = 1, M = 2)
  expect_equal(s, c(-1, 1), tolerance = 1e-12)
})

test_that("make_support: validation errors", {
  expect_error(make_support(-1, M = 5),   "positive")
  expect_error(make_support(1,  M = 1),   "at least 2")
  expect_error(make_support("a"),          "single positive")
})

test_that("default_supports: returns correct structure", {
  set.seed(1)
  n <- 30; K <- 3
  X <- cbind(1, matrix(rnorm(n * (K - 1)), n, K - 1))
  y <- X %*% c(1, 2, -1) + rnorm(n)
  sp <- default_supports(y, X)

  expect_named(sp, c("Z", "V"))
  expect_equal(nrow(sp$Z), K)
  expect_equal(ncol(sp$Z), 5)   # default M_signal = 5
  expect_length(sp$V, 5)         # default M_error = 5
})

test_that("default_supports: error support is symmetric around zero", {
  set.seed(2)
  X <- cbind(1, rnorm(20))
  y <- X %*% c(0.5, 1.5) + rnorm(20)
  sp <- default_supports(y, X)
  expect_equal(sp$V[1], -sp$V[5], tolerance = 1e-12)
  expect_equal(sp$V[3], 0,        tolerance = 1e-12)
})

test_that("default_supports: signal support contains OLS estimate", {
  set.seed(3)
  n <- 50
  X <- cbind(1, rnorm(n))
  beta_true <- c(2, -3)
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  sp <- default_supports(y, X)

  beta_ols <- as.numeric(solve(crossprod(X), crossprod(X, y)))
  for (k in seq_along(beta_ols)) {
    expect_lte(sp$Z[k, 1], beta_ols[k])  # lower bound
    expect_gte(sp$Z[k, 5], beta_ols[k])  # upper bound
  }
})
