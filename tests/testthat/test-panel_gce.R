# Tests for panel_gce() — dual GCE for the one-way error-components panel model
# (Lee & Cheon 2014, eq. 3.12).

# Simulate a balanced panel y_nt = x' beta + mu_n + eps, in the requested layout.
sim_panel <- function(N = 40L, T = 6L, K = 2L, beta = c(1.5, -0.8),
                      seed = 1, layout = "unit") {
  set.seed(seed)
  mu <- rnorm(N)
  if (layout == "unit") id <- rep(seq_len(N), each = T)
  else                  id <- rep(seq_len(N), times = T)
  X <- matrix(rnorm(N * T * K), N * T, K)
  y <- as.vector(X %*% beta) + mu[id] + rnorm(N * T, 0, 0.5)
  list(y = y, X = X, id = id, mu = mu, N = N, T = T, K = K, beta = beta,
       layout = layout)
}

within_beta <- function(d) {
  dm <- function(v) v - ave(v, d$id)
  Xw <- apply(d$X, 2L, dm)
  as.vector(solve(crossprod(Xw), crossprod(Xw, dm(d$y))))
}

wideZ <- function(K) matrix(c(-10, 0, 10), K, 3, byrow = TRUE)
muFF  <- function(N) matrix(c(-6, 0, 6), N, 3, byrow = TRUE)

test_that("converges and satisfies the model FOC", {
  d <- sim_panel()
  fit <- panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N), layout = "unit")

  expect_s3_class(fit, "panel_gce")
  expect_s3_class(fit, "infometrics")
  expect_equal(fit$method, "dual")
  expect_true(fit$converged)
  expect_lt(fit$foc_residual, 1e-3)
})

test_that("recovers the within (fixed-effects) estimate, both layouts", {
  for (lay in c("unit", "period")) {
    d <- sim_panel(layout = lay, seed = if (lay == "unit") 1L else 2L)
    fit <- panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N), layout = lay)
    expect_lt(max(abs(unname(coef(fit)) - within_beta(d))), 0.02)
  }
})

test_that("recovers the individual effects", {
  d <- sim_panel()
  fit <- panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N), layout = "unit")
  expect_gt(cor(fit$mu_hat, d$mu), 0.9)
})

test_that("the prior pulls beta toward its mass; a narrow support shrinks beta", {
  # beta is well-determined here, so (as in markov_ce) the prior barely moves it
  # under a wide support. On a moderate support the prior bites: a prior with
  # mass on the low end of Z pulls beta down, mass on the high end pulls it up.
  d    <- sim_panel()
  Zmod <- matrix(c(-2, 0, 2), d$K, 3, byrow = TRUE)          # moderate support
  p_lo <- matrix(c(.80, .15, .05), d$K, 3, byrow = TRUE)     # mass on the low end
  p_hi <- matrix(c(.05, .15, .80), d$K, 3, byrow = TRUE)     # mass on the high end
  f_lo <- panel_gce(d$y, d$X, Zmod, tt = d$T, FF = muFF(d$N), layout = "unit", p0 = p_lo)
  f_hi <- panel_gce(d$y, d$X, Zmod, tt = d$T, FF = muFF(d$N), layout = "unit", p0 = p_hi)
  expect_true(all(f_hi$b_hat > f_lo$b_hat))                  # prior pulls beta toward its mass

  f_ref    <- panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N), layout = "unit")
  f_narrow <- panel_gce(d$y, d$X, matrix(c(-1, 0, 1), d$K, 3, byrow = TRUE),
                        tt = d$T, FF = muFF(d$N), layout = "unit")
  expect_lt(abs(f_narrow$b_hat[1]), abs(f_ref$b_hat[1]))     # narrow support shrinks toward 0
})

test_that("input validation throws informative errors", {
  d <- sim_panel(N = 20L, T = 4L)
  expect_error(panel_gce(d$y, d$X, wideZ(d$K), tt = 3L, FF = muFF(d$N)),
               "multiple of tt")
  expect_error(panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = matrix(0, d$N + 1, 3)),
               "nrow\\(FF\\)")
  expect_error(panel_gce(d$y, d$X, matrix(0, d$K + 1, 3), tt = d$T, FF = muFF(d$N)),
               "nrow\\(Z\\)")
  expect_error(panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N), nu = 0),
               "inside \\(0, 1\\)")
  expect_error(panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N),
                         p0 = matrix(1 / 2, d$K, 2)), "p0")
})

test_that("draft names, canonical aliases, and S3 methods are consistent", {
  d <- sim_panel(N = 20L, T = 5L)
  fit <- panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N), layout = "unit")

  expect_identical(fit$coefficients, fit$b_hat)
  expect_identical(coef(fit), fit$b_hat)
  expect_identical(fit$objective, fit$value)
  expect_equal(fit$converged, fit$convergence == 0L)
  expect_length(fit$H_signal, d$K)
  expect_length(coef(fit), d$K)
  expect_length(fitted(fit), d$N * d$T)
  expect_length(residuals(fit), d$N * d$T)
  expect_equal(unname(fitted(fit) + residuals(fit)), d$y, tolerance = 1e-10)
  expect_length(fit$mu_hat, d$N)
  expect_output(print(fit), "Panel GCE")
  expect_output(print(summary(fit)), "Coefficients")
})

test_that("standard errors (beta, p, mu) are present, consistent, and finite", {
  d   <- sim_panel()
  fit <- panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N), layout = "unit")

  V <- vcov(fit)
  expect_equal(dim(V), c(d$K, d$K))
  expect_equal(V, t(V))                                        # symmetric
  expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 0))  # PD
  expect_length(fit$se_beta, d$K)
  expect_equal(unname(fit$se_beta), unname(sqrt(diag(V))))

  expect_equal(dim(fit$se_p), c(d$K, ncol(wideZ(d$K))))        # K x M
  expect_true(all(is.finite(fit$se_p)))
  expect_true(all(fit$se_p >= 0))

  expect_length(fit$se_mu, d$N)
  expect_true(all(is.finite(fit$se_mu) & fit$se_mu > 0))
  # mu SE is dominated by the finite-T floor sqrt(sigma2_eps / T)
  expect_true(all(fit$se_mu > sqrt(fit$sigma2_eps / d$T) * 0.9))
})

test_that("vcov_type variants run; disturbance >= within (Sec. 3.3)", {
  d  <- sim_panel()
  fw <- panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N),
                  layout = "unit", vcov_type = "within")
  fc <- panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N),
                  layout = "unit", vcov_type = "cluster")
  fd <- panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N),
                  layout = "unit", vcov_type = "disturbance")
  expect_identical(fw$vcov_type, "within")
  expect_true(all(is.finite(fc$se_beta)) && all(is.finite(fd$se_beta)))
  # the paper's composite-disturbance form is conservative (MC-validated ordering)
  expect_true(all(fd$se_beta >= fw$se_beta - 1e-8))
})

test_that("se_p reproduces the analytic delta formula", {
  d   <- sim_panel()
  Z   <- wideZ(d$K)
  fit <- panel_gce(d$y, d$X, Z, tt = d$T, FF = muFF(d$N), layout = "unit")
  man <- matrix(0, d$K, ncol(Z))
  for (k in seq_len(d$K)) {
    s2       <- sum(fit$p_hat[k, ] * (Z[k, ] - fit$b_hat[k])^2)
    man[k, ] <- abs(fit$p_hat[k, ] * (Z[k, ] - fit$b_hat[k])) / s2 * fit$se_beta[k]
  }
  expect_equal(unname(fit$se_p), man, tolerance = 1e-8)
})

test_that("summary reports a p-value column and print shows SEs", {
  d   <- sim_panel(N = 20L, T = 5L)
  fit <- panel_gce(d$y, d$X, wideZ(d$K), tt = d$T, FF = muFF(d$N), layout = "unit")
  expect_output(print(summary(fit)), "Pr\\(>\\|z\\|\\)")
  expect_output(print(fit), "se")
})
