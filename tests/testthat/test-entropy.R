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

test_that("kl_divergence: D(p, p) = 0", {
  p <- c(0.5, 0.3, 0.2)
  expect_equal(kl_divergence(p, p), 0, tolerance = 1e-12)
})

test_that("kl_divergence: D(p, q) > 0 when p != q", {
  p <- c(0.5, 0.3, 0.2)
  q <- rep(1/3, 3)
  expect_gt(kl_divergence(p, q), 0)
})

test_that("kl_divergence: asymmetry D(p,q) != D(q,p)", {
  p <- c(0.5, 0.3, 0.2)
  q <- rep(1/3, 3)
  expect_false(isTRUE(all.equal(kl_divergence(p, q), kl_divergence(q, p))))
})

test_that("kl_divergence: handles zero in p with 0*log(0/q) = 0", {
  p <- c(0, 0.5, 0.5)
  q <- c(0.2, 0.4, 0.4)
  # Should not produce NaN or Inf
  result <- kl_divergence(p, q)
  expect_true(is.finite(result))
  expect_gte(result, 0)
})

test_that("kl_divergence: error when q = 0 where p > 0", {
  p <- c(0.5, 0.5)
  q <- c(1, 0)
  expect_error(kl_divergence(p, q), "support")
})

test_that("renyi_entropy: converges to Shannon as alpha -> 1", {
  p <- c(0.5, 0.3, 0.2)
  h_shannon <- shannon_entropy(p)
  h_renyi   <- renyi_entropy(p, alpha = 1 - 1e-9)
  expect_equal(h_renyi, h_shannon, tolerance = 1e-5)
})

test_that("renyi_entropy: alpha = 2 gives collision entropy", {
  p <- c(0.5, 0.3, 0.2)
  # H^R_2 = -log(sum(p^2))
  expected <- -log(sum(p^2))
  expect_equal(renyi_entropy(p, alpha = 2), expected, tolerance = 1e-12)
})

test_that("renyi_entropy: uniform gives log(K) for all alpha", {
  K <- 6
  p <- rep(1/K, K)
  for (a in c(0.5, 2, 3)) {
    expect_equal(renyi_entropy(p, alpha = a), log(K), tolerance = 1e-10)
  }
})

test_that("renyi_entropy: rejects alpha = 1", {
  expect_error(renyi_entropy(c(0.5, 0.5), alpha = 1), "alpha = 1")
})

test_that("tsallis_divergence: converges to KL as alpha -> 1", {
  p <- c(0.5, 0.3, 0.2)
  q <- rep(1/3, 3)
  kl <- kl_divergence(p, q)
  ts <- tsallis_divergence(p, q, alpha = 1 - 1e-9)
  expect_equal(ts, kl, tolerance = 1e-5)
})

test_that("cressie_read: Pearson chi-squared at alpha = 1", {
  p <- c(0.5, 0.3, 0.2)
  q <- rep(1/3, 3)
  # CR(alpha=1) = (1/2) * sum(p * ((p/q) - 1))
  # = (1/2) * sum((p - q)^2 / q)  [Pearson chi-sq / 2]
  expected <- (1/2) * sum((p - q)^2 / q)
  expect_equal(cressie_read(p, q, alpha = 1), expected, tolerance = 1e-12)
})

test_that("cressie_read: converges to KL as alpha -> 0", {
  p <- c(0.5, 0.3, 0.2)
  q <- rep(1/3, 3)
  kl <- kl_divergence(p, q)
  cr <- cressie_read(p, q, alpha = 1e-9)
  expect_equal(cr, kl, tolerance = 1e-5)
})

test_that("cressie_read: errors for alpha = 0 and alpha = -1", {
  p <- c(0.5, 0.5); q <- c(0.5, 0.5)
  expect_error(cressie_read(p, q, alpha = 0),  "alpha = 0")
  expect_error(cressie_read(p, q, alpha = -1), "alpha = -1")
})

test_that("normalized_entropy: uniform gives 1", {
  expect_equal(normalized_entropy(rep(1/6, 6)), 1, tolerance = 1e-12)
})

test_that("normalized_entropy: point mass gives 0", {
  expect_equal(normalized_entropy(c(1, 0, 0, 0)), 0, tolerance = 1e-12)
})

test_that("normalized_entropy: values in [0, 1]", {
  set.seed(1)
  for (i in 1:20) {
    x  <- runif(10); p <- x / sum(x)
    sn <- normalized_entropy(p)
    expect_gte(sn, 0 - 1e-12)
    expect_lte(sn, 1 + 1e-12)
  }
})
