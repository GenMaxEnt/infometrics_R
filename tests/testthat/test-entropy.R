test_that("shannon_entropy: uniform distribution gives log(K)", {
  for (K in c(2, 6, 10)) {
    p <- rep(1/K, K)
    expect_equal(shannon_entropy(p), log(K), tolerance = 1e-12)
  }
})

test_that("shannon_entropy: point mass gives 0", {
  p <- c(1, 0, 0, 0)
  expect_equal(shannon_entropy(p), 0, tolerance = 1e-12)
})

test_that("shannon_entropy: binary entropy H(0.5, 0.5) = log(2)", {
  expect_equal(shannon_entropy(c(0.5, 0.5)), log(2), tolerance = 1e-12)
})

test_that("shannon_entropy: base = 2 gives bits", {
  # H(0.5, 0.5) = 1 bit
  expect_equal(shannon_entropy(c(0.5, 0.5), base = 2), 1, tolerance = 1e-12)
})

test_that("shannon_entropy: monotone in uniformity", {
  p_uniform <- rep(1/6, 6)
  p_skewed  <- c(0.5, 0.3, 0.1, 0.05, 0.03, 0.02)
  expect_gt(shannon_entropy(p_uniform), shannon_entropy(p_skewed))
})

test_that("shannon_entropy: input validation", {
  expect_error(shannon_entropy(c(-0.1, 0.6, 0.5)), "non-negative")
  expect_error(shannon_entropy(c(0.5, 0.4)),       "sum to 1")
  expect_error(shannon_entropy("a"),                "numeric")
})
