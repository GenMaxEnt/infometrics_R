# Tests for matrix_ce() — cross-entropy matrix balancing (Golan 2008, Sec. 7.2).

# Build the max-entropy matrix p_ij = exp(lambda_i x_j) / Omega_j (uniform prior)
# from a chosen lambda, so the data y = p %*% x is exactly recoverable.
make_known <- function(lambda, x) {
  A  <- outer(lambda, x)                       # n x m, log p0 = 0 (uniform)
  EA <- exp(sweep(A, 2L, apply(A, 2L, max), "-"))
  sweep(EA, 2L, colSums(EA), "/")
}

test_that("recovers a known maximum-entropy matrix", {
  lambda <- c(0.5, -0.3, 0.1)
  x <- c(2, 5)
  P <- make_known(lambda, x)                   # 3 x 2, columns sum to 1
  y <- as.vector(P %*% x)
  fit <- matrix_ce(y, x)

  expect_s3_class(fit, "matrix_ce")
  expect_s3_class(fit, "infometrics")
  expect_equal(unname(fit$p), P, tolerance = 1e-6)
})

test_that("moment and column-sum constraints are satisfied", {
  set.seed(7)
  Pt <- matrix(runif(12), 4, 3); Pt <- sweep(Pt, 2, colSums(Pt), "/")
  x  <- runif(3, 1, 5)
  y  <- as.vector(Pt %*% x)
  fit <- matrix_ce(y, x)

  expect_lt(max(abs(residuals(fit))), 1e-4)
  expect_lt(max(abs(colSums(coef(fit)) - 1)), 1e-8)
  expect_true(fit$converged)
  expect_equal(fit$method, "dual")
})

test_that("CE with uniform prior equals ME (p0 = NULL)", {
  P <- sweep(matrix(1:6, ncol = 2), 2, colSums(matrix(1:6, ncol = 2)), "/")
  x <- c(3, 2); y <- as.vector(P %*% x)

  fit_me <- matrix_ce(y, x)
  fit_ce <- matrix_ce(y, x, p0 = matrix(1 / 3, 3, 2))
  expect_equal(fit_me$p, fit_ce$p, tolerance = 1e-8)
})

test_that("a non-uniform prior shifts the solution but keeps moments", {
  P <- sweep(matrix(1:6, ncol = 2), 2, colSums(matrix(1:6, ncol = 2)), "/")
  x <- c(3, 2); y <- as.vector(P %*% x)

  fit_unif <- matrix_ce(y, x)
  p0s      <- matrix(c(.6, .3, .1, .1, .3, .6), 3, 2)
  fit_ce   <- matrix_ce(y, x, p0 = p0s)

  expect_gt(max(abs(fit_ce$p - fit_unif$p)), 1e-3)
  expect_lt(max(abs(residuals(fit_ce))), 1e-4)
})

test_that("entropy diagnostics are well-behaved", {
  set.seed(11)
  Pt <- matrix(runif(12), 4, 3); Pt <- sweep(Pt, 2, colSums(Pt), "/")
  x  <- runif(3, 1, 5); y <- as.vector(Pt %*% x)
  fit <- matrix_ce(y, x)

  expect_gte(fit$S, 0); expect_lte(fit$S, 1)
  expect_gte(fit$cross_entropy, 0)            # D(p||uniform) >= 0
  expect_equal(fit$H_signal, fit$entropy)     # canonical alias

  # D(p || p) = 0 when the prior equals the solution
  fit2 <- matrix_ce(y, x, p0 = fit$p)
  expect_lt(fit2$cross_entropy, 1e-6)
})

test_that("input validation throws informative errors", {
  P <- sweep(matrix(1:6, ncol = 2), 2, colSums(matrix(1:6, ncol = 2)), "/")
  x <- c(3, 2); y <- as.vector(P %*% x)

  expect_error(matrix_ce(y, x, p0 = matrix(1 / 2, 2, 2)),
               "dim\\(p0\\) must be")
  expect_error(matrix_ce(y, x, p0 = matrix(c(-1, 1, 1, 1, 1, 1), 3, 2)),
               "strictly positive")
  expect_error(matrix_ce(c("a", "b"), x), "must be numeric")
  expect_error(matrix_ce(c(1, NA, 2), x), "must not contain NA")
})

test_that("a non-summing prior is renormalized with a warning", {
  P <- sweep(matrix(1:6, ncol = 2), 2, colSums(matrix(1:6, ncol = 2)), "/")
  x <- c(3, 2); y <- as.vector(P %*% x)
  expect_warning(
    fit <- matrix_ce(y, x, p0 = matrix(2, 3, 2)),   # columns sum to 6
    "renormaliz"
  )
  expect_equal(colSums(fit$prior), c(1, 1), tolerance = 1e-12)
})

test_that("S3 methods run and return correctly shaped objects", {
  P <- sweep(matrix(1:6, ncol = 2), 2, colSums(matrix(1:6, ncol = 2)), "/")
  x <- c(3, 2); y <- as.vector(P %*% x)
  fit <- matrix_ce(y, x)

  expect_equal(dim(coef(fit)), c(3L, 2L))     # n x m matrix
  expect_length(fitted(fit), 3)
  expect_length(residuals(fit), 3)
  expect_output(print(fit), "matrix balancing")
  s <- summary(fit)
  expect_s3_class(s, "summary.matrix_ce")
  expect_output(print(s), "Moment fit")
})

test_that("normalized entropy is prior-relative and matches the uniform max", {
  P <- sweep(matrix(1:6, ncol = 2), 2, colSums(matrix(1:6, ncol = 2)), "/")
  x <- c(3, 2); y <- as.vector(P %*% x)
  fit <- matrix_ce(y, x)                         # uniform p0 (ME)

  # For a uniform prior, H(p0) = m * log(n), so S = H(p)/(m log n) (unchanged).
  expect_equal(fit$entropy_p0, ncol(fit$p) * log(nrow(fit$p)), tolerance = 1e-10)
  expect_equal(fit$S, fit$entropy / (ncol(fit$p) * log(nrow(fit$p))),
               tolerance = 1e-10)
  expect_gte(fit$S, 0); expect_lte(fit$S, 1)
})

test_that("summary reports the Golan 7.5 information measures", {
  P <- sweep(matrix(1:6, ncol = 2), 2, colSums(matrix(1:6, ncol = 2)), "/")
  x <- c(3, 2); y <- as.vector(P %*% x)
  s <- summary(matrix_ce(y, x))

  expect_true(all(c("S", "I", "pseudo_R2", "S_col", "W", "df", "p_value")
                  %in% names(s)))
  expect_equal(s$I, 1 - s$S)
  expect_equal(s$pseudo_R2, 1 - s$S)
  expect_length(s$S_col, 2L)                     # one per column
  expect_equal(as.integer(s$df), 2L)             # n - 1 = 3 - 1
  expect_true(is.finite(s$p_value) && s$p_value >= 0 && s$p_value <= 1)
  expect_gte(s$W, 0)                             # data moves p off the uniform prior
  expect_output(print(s), "Information measures")
})

test_that("entropy-ratio test is ~0 when the data is uninformative", {
  # y_i = (1/n) sum_j x_j  =>  p = p0 (uniform) exactly  =>  W = 0, p ~ 1
  x <- c(1, 4, 9, 16); y <- rep(mean(x), 3)
  suppressWarnings(s <- summary(matrix_ce(y, x)))
  expect_lt(s$W, 1e-6)
  expect_gt(s$p_value, 0.99)
})

test_that("uniform-consistent data recovers the uniform matrix", {
  # y_i = (1/n) sum_j x_j for all i  =>  p = 1/n exactly. The optimum is
  # lambda = 0, where optim can report a benign convergence false-negative.
  x <- c(1, 4, 9, 16)
  y <- rep(mean(x), 3)
  suppressWarnings(fit <- matrix_ce(y, x))
  expect_lt(max(abs(fit$p - 1 / 3)), 1e-6)
})
