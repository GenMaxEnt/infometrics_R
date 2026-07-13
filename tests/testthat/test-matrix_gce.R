# Tests for matrix_gce() — GCE matrix balancing with stochastic moments
# (Golan 2008, Sec. 7.4). Inputs are normalized to [0, 1].

make_gce_data <- function() {
  set.seed(1)
  P <- matrix(runif(6), ncol = 2); P <- sweep(P, 2, colSums(P), "/")
  x <- c(0.6, 0.3)
  list(P = P, x = x, y = as.vector(P %*% x))
}

test_that("converges and satisfies the stochastic moment condition", {
  d <- make_gce_data()
  fit <- matrix_gce(d$y, d$x, nu = 0.5)

  expect_s3_class(fit, "matrix_gce")
  expect_s3_class(fit, "infometrics")
  expect_true(fit$converged)
  expect_equal(fit$method, "dual")
  # y - Px = e at the optimum
  expect_lt(max(abs(residuals(fit) - fit$e)), 1e-6)
  expect_lt(max(abs(colSums(coef(fit)) - 1)), 1e-8)   # signal columns sum to 1
  expect_lt(max(abs(rowSums(fit$w) - 1)), 1e-8)       # noise rows sum to 1
})

test_that("as the noise support shrinks, the signal returns to matrix_ce", {
  d <- make_gce_data()
  fit_gce <- matrix_gce(d$y, d$x, v = c(-1e-4, 0, 1e-4), nu = 0.5)
  fit_ce  <- matrix_ce(d$y, d$x)
  expect_equal(unname(coef(fit_gce)), unname(coef(fit_ce)), tolerance = 1e-4)
})

test_that("GME (uniform priors) differs from GCE (non-uniform signal prior)", {
  d <- make_gce_data()
  fit_gme <- matrix_gce(d$y, d$x, nu = 0.5)
  p0s     <- matrix(c(.6, .3, .1, .1, .3, .6), 3, 2)
  fit_gce <- matrix_gce(d$y, d$x, p0 = p0s, nu = 0.5)

  expect_gt(max(abs(fit_gce$p - fit_gme$p)), 1e-3)
  expect_lt(max(abs(residuals(fit_gce) - fit_gce$e)), 1e-6)
})

test_that("entropy diagnostics are well-behaved", {
  d <- make_gce_data()
  fit <- matrix_gce(d$y, d$x, nu = 0.5)

  expect_gte(fit$S_p, 0); expect_lte(fit$S_p, 1)
  expect_gte(fit$S_w, 0); expect_lte(fit$S_w, 1)
  expect_equal(fit$H_signal, fit$entropy_p)     # canonical alias
  expect_equal(fit$S, fit$S_p)                  # canonical alias
})

test_that("input validation throws informative errors", {
  d <- make_gce_data()

  expect_error(matrix_gce(d$y, d$x, nu = 0),   "strictly between 0 and 1")
  expect_error(matrix_gce(d$y, d$x, nu = 1),   "strictly between 0 and 1")
  expect_error(matrix_gce(d$y, d$x, nu = 1.5), "strictly between 0 and 1")
  expect_error(matrix_gce(d$y, d$x, p0 = matrix(1 / 2, 2, 2)), "dim\\(p0\\) must be")
  expect_error(matrix_gce(d$y, d$x, v = 1L),   "must be >= 2")
  expect_error(
    matrix_gce(d$y, d$x, v = c(-1, 0, 1), w0 = matrix(1 / 2, 3, 2)),
    "ncol\\(w0\\).*must equal length\\(v\\)"
  )
  expect_error(
    matrix_gce(d$y, d$x, v = c(-1, 0, 1), w0 = matrix(1 / 3, 2, 3)),
    "nrow\\(w0\\).*must equal length\\(y\\)"
  )
})

test_that("off-[0,1] data triggers the normalization warning", {
  expect_warning(matrix_gce(c(2, 5, 8), c(3, 4)), "normalized to \\[0, 1\\]")
})

test_that("S3 methods run and return correctly shaped objects", {
  d <- make_gce_data()
  fit <- matrix_gce(d$y, d$x, nu = 0.5)

  expect_equal(dim(coef(fit)), c(3L, 2L))       # n x m signal matrix
  expect_length(fitted(fit), 3)
  expect_length(residuals(fit), 3)
  expect_output(print(fit), "GCE matrix balancing")

  s <- summary(fit)
  expect_s3_class(s, "summary.matrix_gce")
  expect_true(all(c("W", "df", "p_value", "S_col") %in% names(s)))
  expect_equal(as.integer(s$df), 2L)            # n - 1
  expect_true(is.finite(s$p_value) && s$p_value >= 0 && s$p_value <= 1)
  expect_output(print(s), "Information measures")
})
