# Tests for linreg() — GCE linear regression via the formula interface.

make_reg_data <- function() {
  set.seed(1)
  d <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  d$y <- 2 - 1.5 * d$x1 + 0.8 * d$x2 + rnorm(80, sd = 0.7)
  d
}

test_that("converges and returns a well-formed object", {
  fit <- linreg(y ~ x1 + x2, data = make_reg_data(),
                Z = seq(-20, 20, length.out = 5))

  expect_s3_class(fit, "linreg")
  expect_s3_class(fit, "infometrics")
  expect_true(fit$converged)
  expect_equal(fit$method, "dual")
  expect_length(coef(fit), 3)                  # intercept + x1 + x2
  expect_named(coef(fit), c("(Intercept)", "x1", "x2"))
})

test_that("GCE approaches OLS as the supports widen", {
  d   <- make_reg_data()
  ols <- coef(lm(y ~ x1 + x2, data = d))

  fit_wide  <- linreg(y ~ x1 + x2, data = d,
                      Z = seq(-200, 200, length.out = 5),
                      v = seq(-200, 200, length.out = 3))
  fit_tight <- linreg(y ~ x1 + x2, data = d,
                      Z = seq(-3, 3, length.out = 5))

  # Wide supports => close to OLS; and closer than tight (which shrinks to 0).
  expect_lt(max(abs(coef(fit_wide) - ols)), 0.05)
  expect_lt(max(abs(coef(fit_wide) - ols)),
            max(abs(coef(fit_tight) - ols)))
})

test_that("tight supports shrink coefficients toward the support center (0)", {
  d <- make_reg_data()
  fit_tight <- linreg(y ~ x1 + x2, data = d, Z = seq(-3, 3, length.out = 5))
  fit_wide  <- linreg(y ~ x1 + x2, data = d, Z = seq(-200, 200, length.out = 5))
  # The slope magnitudes are pulled in under tight supports.
  expect_lt(abs(coef(fit_tight)["x1"]), abs(coef(fit_wide)["x1"]))
})

test_that("moment condition holds: residuals equal the estimated noise", {
  fit <- linreg(y ~ x1 + x2, data = make_reg_data(),
                Z = seq(-20, 20, length.out = 5))
  # grad = Xbeta + e - y = 0  =>  y - Xbeta = e.
  expect_lt(max(abs(residuals(fit) - fit$e)), 1e-4)
})

test_that("GME (uniform p0) differs from GCE (non-uniform signal prior)", {
  d  <- make_reg_data()
  Zc <- seq(-20, 20, length.out = 5)
  fit_gme <- linreg(y ~ x1 + x2, data = d, Z = Zc)              # uniform p0
  p0_skew <- matrix(c(.4, .3, .15, .1, .05), nrow = 3, ncol = 5, byrow = TRUE)
  fit_gce <- linreg(y ~ x1 + x2, data = d, Z = Zc, p0 = p0_skew)
  expect_false(isTRUE(all.equal(unname(coef(fit_gme)),
                                unname(coef(fit_gce)), tolerance = 1e-3)))
})

test_that("a shared support vector and the equivalent matrix agree", {
  d  <- make_reg_data()
  zv <- seq(-15, 15, length.out = 5)
  zm <- matrix(zv, nrow = 3, ncol = 5, byrow = TRUE)
  fit_v <- linreg(y ~ x1 + x2, data = d, Z = zv)
  fit_m <- linreg(y ~ x1 + x2, data = d, Z = zm)
  expect_equal(coef(fit_v), coef(fit_m), tolerance = 1e-8)
})

test_that("predict returns fitted values and handles newdata", {
  d   <- make_reg_data()
  fit <- linreg(y ~ x1 + x2, data = d, Z = seq(-20, 20, length.out = 5))

  expect_equal(predict(fit), fitted(fit), tolerance = 1e-12)
  nd <- data.frame(x1 = c(0, 1, -1), x2 = c(0, -1, 1))
  pr <- predict(fit, newdata = nd)
  expect_length(pr, 3)
  # Manual check against coefficients
  b <- coef(fit)
  expect_equal(pr, as.vector(b[1] + b[2] * nd$x1 + b[3] * nd$x2),
               tolerance = 1e-10)
})

test_that("vcov is symmetric positive definite with the right dimensions", {
  fit <- linreg(y ~ x1 + x2, data = make_reg_data(),
                Z = seq(-20, 20, length.out = 5))
  V <- vcov(fit)
  expect_equal(dim(V), c(3L, 3L))
  expect_equal(V, t(V), tolerance = 1e-12)
  expect_gt(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_true(all(is.finite(sqrt(diag(V)))))
})

test_that("analytic dual Hessian is symmetric positive definite", {
  fit <- linreg(y ~ x1 + x2, data = make_reg_data(),
                Z = seq(-20, 20, length.out = 5))
  H <- fit$hessian
  expect_equal(H, t(H), tolerance = 1e-8)
  expect_gt(min(eigen(H, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("entropy S in [0, 1] and r.squared in a sensible range", {
  fit <- linreg(y ~ x1 + x2, data = make_reg_data(),
                Z = seq(-20, 20, length.out = 5))
  expect_gte(fit$S, 0)
  expect_lte(fit$S, 1)
  expect_gt(fit$r.squared, 0.5)
  expect_lte(fit$r.squared, 1)
  expect_length(fit$H_signal, 3)               # per-coefficient signal entropy
})

test_that("input validation throws informative errors", {
  d <- make_reg_data()

  expect_error(linreg(y ~ x1 + x2, data = d, nu = 0),
               "strictly between 0 and 1")
  expect_error(linreg(y ~ x1 + x2, data = d, nu = 1),
               "strictly between 0 and 1")
  # Z with wrong number of rows (K = 3 here)
  expect_error(
    linreg(y ~ x1 + x2, data = d, Z = matrix(seq(-5, 5, length.out = 10), 2, 5)),
    "must equal the number of model-matrix columns"
  )
  # p0 dimension mismatch
  expect_error(
    linreg(y ~ x1 + x2, data = d, Z = seq(-5, 5, length.out = 5),
           p0 = matrix(1 / 4, 3, 4)),
    "must match dim\\(Z\\)"
  )
  # scalar v < 2
  expect_error(linreg(y ~ x1 + x2, data = d, v = 1L), "must be >= 2")
  # w0 with wrong number of columns
  expect_error(
    linreg(y ~ x1 + x2, data = d, v = c(-1, 0, 1), w0 = matrix(1 / 2, 80, 2)),
    "must equal length\\(v\\)"
  )
})

test_that("S3 methods run and return correctly shaped objects", {
  fit <- linreg(y ~ x1 + x2, data = make_reg_data(),
                Z = seq(-20, 20, length.out = 5))

  expect_length(coef(fit), 3)
  expect_length(fitted(fit), 80)
  expect_length(residuals(fit), 80)
  expect_equal(dim(vcov(fit)), c(3L, 3L))
  expect_output(print(fit), "Coefficients")
  s <- summary(fit)
  expect_s3_class(s, "summary.linreg")
  expect_output(print(s), "Coefficients")
  expect_equal(dim(s$coefficients), c(3L, 4L))   # Est, SE, ER, Pr(>Chi)
  expect_equal(colnames(s$coefficients),
               c("Estimate", "Std. Error", "ER", "Pr(>Chi)"))
})

test_that("summary reports a chi-squared entropy-ratio test with valid stats", {
  fit <- linreg(y ~ x1 + x2, data = make_reg_data(),
                Z = seq(-20, 20, length.out = 5))
  s  <- summary(fit)
  ER <- s$coefficients[, "ER"]
  pv <- s$coefficients[, "Pr(>Chi)"]

  expect_true(all(ER >= 0))                       # clamped, non-negative
  expect_true(all(pv >= 0 & pv <= 1))             # valid p-values
  expect_lt(pv["x1"], 0.05)                       # strong effect rejects H0
  # ER row exists for every coefficient, including the intercept
  expect_length(ER, 3L)
  expect_false(is.na(ER["(Intercept)"]))
})

test_that("ER test recovers true zero vs non-zero coefficients (Golan 6.6 style)", {
  set.seed(11)
  n  <- 300
  x2 <- runif(n, 0, 20)
  x3 <- runif(n, 0, 20)
  x4 <- x3 + rnorm(n, 0, 0.2)                     # near-collinear with x3 (true beta4 = 0)
  y  <- 1 - 2 * x2 + 3 * x3 + 0 * x4 + rnorm(n, 0, sqrt(2))
  d  <- data.frame(y = y, x2 = x2, x3 = x3, x4 = x4)
  Zs <- matrix(rep(c(-100, -50, 0, 50, 100), each = 4), nrow = 4)   # 4 x 5

  s  <- summary(linreg(y ~ x2 + x3 + x4, data = d, Z = Zs))
  pv <- s$coefficients[, "Pr(>Chi)"]

  # beta2 is a clean (non-collinear) strong effect -> the ER test rejects.
  expect_lt(pv["x2"], 0.05)
  # beta4 is truly 0 -> not rejected.
  expect_gt(pv["x4"], 0.10)
  # x3 and x4 are near-collinear (x4 = x3 + small noise), so an INDIVIDUAL ER
  # test for beta3 = 0 need not reject -- x4 absorbs the dropped x3 effect.
  # This is exactly the book's collinearity illustration; the joint/overall
  # test is what detects the signal. So we assert the overall test, not an
  # individual x3 rejection.
  # Overall ER test: H0 all coefficients (intercept + 3 slopes) = 0 -> df = K = 4.
  expect_lt(s$er_overall[["p.value"]], 0.05)      # coefficients jointly non-zero -> reject
  expect_equal(as.integer(s$er_overall[["df"]]), 4L)
  expect_output(print(s), "Overall ER")
})

test_that("overall ER test stays finite when the restricted model degenerates", {
  # Strong level, low variance, (near) noise-free: y has a large mean relative
  # to sd(y), so the zero-centred noise support cannot span y once the overall
  # test forces the intercept to 0. The restricted dual is then unbounded and
  # lambda diverges, underflowing the noise softmax to exact zeros. H* must
  # apply the 0*log0 = 0 convention so the statistic stays finite (and the line
  # is reported) rather than collapsing to NaN.
  set.seed(125)
  XX <- cbind(1, matrix(runif(20), ncol = 2))
  yy <- as.vector(XX %*% c(3, 2, 5))
  d  <- data.frame(x2 = XX[, 2], x3 = XX[, 3], y = yy)
  ZZ <- rbind(c(-10, 0, 10), c(-10, 0, 10), c(-10, 0, 10))

  s  <- summary(linreg(y ~ x2 + x3, data = d, Z = ZZ))
  eo <- s$er_overall
  expect_true(is.finite(eo[["ER"]]))             # not NaN
  expect_true(is.finite(eo[["p.value"]]))
  expect_gte(eo[["ER"]], 0)
  expect_equal(as.integer(eo[["df"]]), 3L)       # all coefficients (intercept + 2)
  expect_lt(eo[["p.value"]], 0.05)               # coefficients clearly non-zero -> reject
  expect_output(print(s), "Overall ER")          # line is actually printed
})
