# Tests for multinomial_gce() — nu-weighted GCE multinomial (GJP 1996).

make_mnl <- function(n = 400L, seed = 123) {
  set.seed(seed)
  J <- 3L; x <- rnorm(n); X <- cbind(1, x)
  bt <- cbind(c(0, 0), c(1, 2), c(-0.5, -1))
  P  <- exp(X %*% bt); P <- P / rowSums(P)
  yc <- apply(P, 1, function(p) sample(J, 1, prob = p))
  Y  <- matrix(0, n, J); Y[cbind(seq_len(n), yc)] <- 1
  list(X = X, Y = Y, P = P, yc = yc, x = x, n = n, J = J)
}

test_that("converges, row-stochastic, gradient ~ 0 at optimum", {
  d <- make_mnl()
  fit <- multinomial_gce(d$Y, d$X, which_alternative = 1L, nu = 0.5)

  expect_s3_class(fit, "multinomial_gce")
  expect_s3_class(fit, "infometrics")
  expect_equal(fit$method, "dual")
  expect_true(fit$converged)
  expect_lt(max(abs(rowSums(fit$p_hat) - 1)), 1e-10)
  gr <- infometrics:::.multinomial_gce_grad(
    fit$lambda[, -1], 1L, ncol(d$X), d$n, d$J, d$Y, d$X, fit$p0, fit$w0, fit$v, 0.5)
  expect_lt(max(abs(gr)), 1e-3)
})

test_that("reproduces gme_mnl exactly at nu=0.5 with a matched support", {
  d <- make_mnl()
  vN <- seq(-1 / sqrt(d$n), 1 / sqrt(d$n), length.out = 3)
  f <- multinomial_gce(d$Y, d$X, 1L, v = vN, nu = 0.5)
  g <- gme_mnl(factor(d$yc) ~ x, data = data.frame(yc = d$yc, x = d$x), v = vN)
  expect_lt(max(abs(f$p_hat - g$p_hat)), 1e-5)         # correctness anchor (BFGS tol)
})

test_that("recovers the data-generating probabilities", {
  d <- make_mnl()
  vN <- seq(-1 / sqrt(d$n), 1 / sqrt(d$n), length.out = 3)
  f <- multinomial_gce(d$Y, d$X, 1L, v = vN, nu = 0.5)
  expect_lt(max(abs(f$p_hat - d$P)), 0.1)
})

test_that("nu and support width pull toward the uniform prior", {
  d <- make_mnl()
  meanH <- function(f) mean(f$H_signal)
  f_lo <- multinomial_gce(d$Y, d$X, 1L, v = c(-1, 0, 1), nu = 0.2)
  f_hi <- multinomial_gce(d$Y, d$X, 1L, v = c(-1, 0, 1), nu = 0.8)
  expect_gt(meanH(f_hi), meanH(f_lo))                  # larger nu -> more uniform
})

test_that("a non-uniform signal prior shifts the fit (GCE vs GME)", {
  d <- make_mnl()
  vN <- seq(-1 / sqrt(d$n), 1 / sqrt(d$n), length.out = 3)
  f_gme <- multinomial_gce(d$Y, d$X, 1L, v = vN, nu = 0.5)
  p0s   <- matrix(c(.6, .2, .2), d$n, 3, byrow = TRUE)
  f_gce <- multinomial_gce(d$Y, d$X, 1L, v = vN, nu = 0.5, p0 = p0s)
  expect_gt(max(abs(f_gce$p_hat - f_gme$p_hat)), 1e-3)
})

test_that("input validation throws informative errors", {
  d <- make_mnl(n = 60L)
  badY <- d$Y; badY[1, 1] <- 0.5
  expect_error(multinomial_gce(badY, d$X), "sum to 1")
  expect_error(multinomial_gce(d$Y, d$X[-1, ]), "must equal nrow")
  expect_error(multinomial_gce(d$Y, d$X, nu = 0), "in \\(0, 1\\)")
  expect_error(multinomial_gce(d$Y, d$X, nu = 1), "in \\(0, 1\\)")
  expect_error(multinomial_gce(d$Y, d$X, v = c(-1, 0, 2)), "within \\[-1, 1\\]")
  expect_error(multinomial_gce(d$Y, d$X, p0 = matrix(1 / 2, 2, 2)), "dim\\(p0\\)")
  expect_error(multinomial_gce(d$Y, d$X, which_alternative = 5L),
               "between 1 and")
})

test_that("which_alternative shifts the fit only slightly (GCE noise asymmetry)", {
  d <- make_mnl()
  vN <- seq(-1 / sqrt(d$n), 1 / sqrt(d$n), length.out = 3)
  p1 <- multinomial_gce(d$Y, d$X, 1L, v = vN, nu = 0.5)$p_hat
  p2 <- multinomial_gce(d$Y, d$X, 2L, v = vN, nu = 0.5)$p_hat
  p3 <- multinomial_gce(d$Y, d$X, 3L, v = vN, nu = 0.5)$p_hat
  expect_lt(max(abs(p1 - p2)), 0.05)
  expect_lt(max(abs(p1 - p3)), 0.05)
})

test_that("draft names and canonical aliases are present and consistent", {
  d <- make_mnl(n = 120L)
  fit <- multinomial_gce(d$Y, d$X, 1L, nu = 0.5)
  expect_identical(fit$p_hat, fit$p)
  expect_identical(fit$w_hat, fit$w)
  expect_identical(fit$lambda_hat, fit$lambda)
  expect_identical(fit$e, fit$ee)
  expect_identical(fit$objective, fit$value)
  expect_equal(fit$converged, fit$convergence == 0L)
  expect_equal(fit$method, "dual")
  expect_length(fit$H_signal, d$n)
  expect_identical(coef(fit), fit$lambda)              # B12: coef -> lambda
})

test_that("S3 methods run and return correctly shaped objects", {
  d <- make_mnl(n = 150L)
  fit <- multinomial_gce(d$Y, d$X, 1L, nu = 0.5)
  expect_equal(dim(coef(fit)), c(2L, 3L))              # K x J
  expect_equal(dim(fitted(fit)), c(150L, 3L))
  expect_equal(dim(residuals(fit)), c(150L, 3L))
  expect_output(print(fit), "Multinomial")
  expect_s3_class(summary(fit), "multinomial_gce")     # summary returns the object invisibly
  expect_output(print(summary(fit)), "Alternative")
})

test_that("margins matches a numerical derivative", {
  d <- make_mnl(n = 100L)
  fit <- multinomial_gce(d$Y, d$X, 1L, nu = 0.5)
  ME <- margins(fit, average = FALSE)                  # N x J x K (obs, alt, predictor)

  # K x J average effects = transpose of the obs-mean (J x K) of the per-obs array
  expect_equal(margins(fit), t(apply(ME, c(2, 3), mean)),
               tolerance = 1e-10, ignore_attr = TRUE)

  # numerical derivative w.r.t. column 2 (x) at the fitted lambda / p0
  p_at <- function(Xm) {
    a <- log(fit$p0) + (Xm %*% fit$lambda) / fit$nu
    ea <- exp(a - apply(a, 1L, max)); ea / rowSums(ea)
  }
  h <- 1e-6; Xp <- d$X; Xp[, 2] <- Xp[, 2] + h
  num <- (p_at(Xp) - p_at(d$X)) / h                    # N x J: dp_ij/dx_i2
  expect_lt(max(abs(num - ME[, , 2])), 1e-5)
})

# ---- Fano error bounds (Golan 2008, sec 7.5) --------------------------------

test_that("fano_bounds returns per-observation bounds that hold", {
  d   <- make_mnl()
  fit <- multinomial_gce(d$Y, d$X, 1L, nu = 0.5)
  fb  <- fano_bounds(fit)

  expect_s3_class(fb, "data.frame")
  expect_equal(nrow(fb), d$n)
  expect_named(fb, c("p_max", "pe", "H", "S", "pe_lower"))
  # modal error and normalized entropy are exactly as defined
  expect_equal(fb$pe, 1 - apply(fit$p_hat, 1, max), tolerance = 1e-12)
  expect_true(all(fb$S >= 0 & fb$S <= 1 + 1e-12))
  # weak Fano lower bound holds for every observation
  expect_true(all(fb$pe >= fb$pe_lower - 1e-9))
  # strong Fano inequality holds: H(pe) + pe log(J-1) >= H(p_i)
  Hbin <- function(pe) ifelse(pe > 0 & pe < 1,
                              -pe * log(pe) - (1 - pe) * log(1 - pe), 0)
  expect_true(all(Hbin(fb$pe) + fb$pe * log(d$J - 1) >= fb$H - 1e-9))
  ov <- attr(fb, "overall")
  expect_named(ov, c("mean_pe", "mean_pe_lower", "S_system"))
  expect_gte(ov["mean_pe"], ov["mean_pe_lower"])
})

# ---- marginal-effect standard errors ----------------------------------------

# numerical Jacobian of vec(average ME) in the free multipliers (test helper)
.me_jacobian <- function(fit) {
  wa <- fit$which_alternative; K <- fit$K; J <- fit$J
  lf <- as.vector(fit$lambda[, -wa, drop = FALSE])
  hh <- 1e-5; nf <- length(lf); Jac <- matrix(0, K * J, nf)
  for (mm in seq_len(nf)) { e <- numeric(nf); e[mm] <- hh
    Jac[, mm] <- (infometrics:::.multinomial_gce_me_vec(lf + e, fit$X, fit$p0, fit$nu, wa, K, J) -
                  infometrics:::.multinomial_gce_me_vec(lf - e, fit$X, fit$p0, fit$nu, wa, K, J)) / (2 * hh)
  }
  Jac
}

test_that("margins(se = TRUE) defaults to the robust sandwich", {
  d   <- make_mnl()
  fit <- multinomial_gce(d$Y, d$X, 1L, nu = 0.5)
  m   <- margins(fit, se = TRUE)                       # default se_method = "sandwich"

  expect_s3_class(m, "margins_gce")
  expect_equal(m$se_method, "sandwich")
  expect_equal(m$estimate, margins(fit), tolerance = 1e-10)   # matches the point estimate
  expect_equal(dim(m$se), c(fit$K, fit$J))
  expect_true(all(is.finite(m$se) & m$se >= 0))
  expect_equal(m$z, m$estimate / m$se, tolerance = 1e-12)
  expect_output(print(m), "marginal effects")

  # reproduce the sandwich SE independently
  Jac    <- .me_jacobian(fit)
  se_ref <- matrix(sqrt(pmax(diag(Jac %*% infometrics:::.multinomial_gce_sandwich(fit) %*% t(Jac)), 0)),
                   fit$K, fit$J)
  expect_equal(unname(m$se), unname(se_ref), tolerance = 1e-6)
})

test_that("se_method = 'delta' reproduces the naive Hessian-inverse SE", {
  d   <- make_mnl()
  fit <- multinomial_gce(d$Y, d$X, 1L, nu = 0.5)
  md  <- margins(fit, se = TRUE, se_method = "delta")
  Jac <- .me_jacobian(fit)
  se_ref <- matrix(sqrt(pmax(diag(Jac %*% fit$vcov %*% t(Jac)), 0)), fit$K, fit$J)
  expect_equal(unname(md$se), unname(se_ref), tolerance = 1e-6)
  # the naive delta is conservative: strictly larger than the sandwich
  ms <- margins(fit, se = TRUE)                        # sandwich
  expect_true(all(md$se > ms$se))
})

test_that("margins(se = FALSE) is unchanged; se needs average = TRUE", {
  d   <- make_mnl(n = 150L)
  fit <- multinomial_gce(d$Y, d$X, 1L, nu = 0.5)
  expect_true(is.matrix(margins(fit)))                 # backward compatible
  expect_equal(dim(margins(fit, average = FALSE)), c(fit$N, fit$J, fit$K))
  expect_error(margins(fit, se = TRUE, average = FALSE), "average = TRUE")
})

test_that("sandwich margins SEs roughly match a bootstrap", {
  d   <- make_mnl(n = 150L)
  fit <- multinomial_gce(d$Y, d$X, 1L, nu = 0.5)
  ms  <- margins(fit, se = TRUE)                       # sandwich
  set.seed(11)
  mb  <- margins(fit, se = TRUE, se_method = "bootstrap", B = 40L)
  expect_equal(mb$se_method, "bootstrap")
  expect_equal(dim(mb$se), c(fit$K, fit$J))
  expect_true(all(is.finite(mb$se) & mb$se >= 0))
  expect_true(max(ms$se / mb$se) < 2 && min(ms$se / mb$se) > 0.5)  # same order
})
