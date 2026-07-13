# Tests for inverse_ce() — pure-moment ME/CE via the formula interface,
# including information-matrix SEs (Var(lambda) = I^-1) and Fano error bounds.

test_that("recovers known probabilities from exact moments", {
  # Under-determined, full-rank problem (K = 5 > T = 2). p_true is the
  # max-entropy distribution for its own moments, so it is recovered exactly.
  X <- rbind(c(1, 2, 3, 4, 5),
             c(1, 4, 9, 16, 25))
  p_true <- exp(as.vector(crossprod(X, c(0.2, -0.05))))
  p_true <- p_true / sum(p_true)
  dat <- data.frame(y  = as.vector(X %*% p_true),
                    s1 = X[, 1], s2 = X[, 2], s3 = X[, 3],
                    s4 = X[, 4], s5 = X[, 5])
  fit <- inverse_ce(y ~ s1 + s2 + s3 + s4 + s5 - 1, data = dat)

  expect_s3_class(fit, "inverse_ce")
  expect_s3_class(fit, "infometrics")
  expect_equal(unname(coef(fit)), p_true, tolerance = 1e-6)
  expect_lt(max(abs(residuals(fit))), 1e-6)
})

test_that("output is de-duplicated (canonical fields kept, aliases dropped)", {
  X <- rbind(c(1, 2, 3, 4, 5),
             c(1, 4, 9, 16, 25))
  p_true <- exp(as.vector(crossprod(X, c(0.2, -0.05))))
  p_true <- p_true / sum(p_true)
  dat <- data.frame(y  = as.vector(X %*% p_true),
                    s1 = X[, 1], s2 = X[, 2], s3 = X[, 3],
                    s4 = X[, 4], s5 = X[, 5])
  fit <- inverse_ce(y ~ s1 + s2 + s3 + s4 + s5 - 1, data = dat)
  nm  <- names(fit)
  # dropped duplicate/redundant fields
  expect_false(any(c("p", "lambda", "value", "prior", "q", "converged") %in% nm))
  # kept canonical fields (prior stored as p0)
  expect_true(all(c("p_hat", "lambda_hat", "residuals", "objective", "p0",
                    "convergence") %in% nm))
  # accessors read the kept fields
  expect_equal(unname(coef(fit)), unname(fit$p_hat))
})

test_that("moment constraints are satisfied at the optimum", {
  X <- matrix(c(1, 2, 3, 4, 5,
                2, 1, 4, 3, 5), nrow = 2, byrow = TRUE)
  p_true <- c(.30, .25, .20, .15, .10)
  dat <- data.frame(y = as.vector(X %*% p_true),
                    x1 = X[, 1], x2 = X[, 2], x3 = X[, 3],
                    x4 = X[, 4], x5 = X[, 5])
  fit <- inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat)

  expect_lt(max(abs(residuals(fit))), 1e-4)
  expect_equal(fit$convergence, 0L)
  expect_equal(sum(coef(fit)), 1, tolerance = 1e-10)
})

test_that("CE with uniform prior equals ME (default)", {
  X <- matrix(c(1, 2, 3, 4, 5,
                5, 4, 3, 2, 1), nrow = 2, byrow = TRUE)
  p_true <- c(.30, .25, .20, .15, .10)
  dat <- data.frame(y = as.vector(X %*% p_true),
                    x1 = X[, 1], x2 = X[, 2], x3 = X[, 3],
                    x4 = X[, 4], x5 = X[, 5])

  fit_me <- inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat)
  fit_ce <- inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat,
                       p0 = rep(1 / 5, 5))

  expect_equal(coef(fit_me), coef(fit_ce), tolerance = 1e-8)
})

test_that("non-uniform prior pulls the solution toward the prior", {
  X <- matrix(c(1, 2, 3, 4, 5), nrow = 1)   # single mean constraint, K = 5
  y <- 3.0                                   # = fair-die mean of 1:5
  dat <- data.frame(y = y, x1 = 1, x2 = 2, x3 = 3, x4 = 4, x5 = 5)

  # Uniform prior + mean 3 => uniform p (symmetric)
  fit_unif <- inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat)
  expect_equal(unname(coef(fit_unif)), rep(1 / 5, 5), tolerance = 1e-6)

  # Skewed prior shifts the CE solution away from uniform
  p0_skew <- c(.40, .25, .15, .12, .08)
  fit_ce <- inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat,
                       p0 = p0_skew)
  expect_false(isTRUE(all.equal(unname(coef(fit_ce)), rep(1 / 5, 5),
                                tolerance = 1e-3)))
  expect_lt(max(abs(residuals(fit_ce))), 1e-4)
})

test_that("analytic Hessian is symmetric and positive semidefinite", {
  X <- matrix(c(1, 2, 3, 4, 5,
                5, 4, 3, 2, 1), nrow = 2, byrow = TRUE)
  p_true <- c(.30, .25, .20, .15, .10)
  dat <- data.frame(y = as.vector(X %*% p_true),
                    x1 = X[, 1], x2 = X[, 2], x3 = X[, 3],
                    x4 = X[, 4], x5 = X[, 5])
  fit <- inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat)

  H <- fit$hessian
  expect_equal(H, t(H), tolerance = 1e-10)
  expect_gte(min(eigen(H, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
})

test_that("analytic Hessian equals the finite-difference information matrix", {
  X <- rbind(c(1, 2, 3, 4, 5),
             c(1, 4, 9, 16, 25))
  p_true <- exp(as.vector(crossprod(X, c(0.2, -0.05))))
  p_true <- p_true / sum(p_true)
  dat <- data.frame(y  = as.vector(X %*% p_true),
                    s1 = X[, 1], s2 = X[, 2], s3 = X[, 3],
                    s4 = X[, 4], s5 = X[, 5])
  fit <- inverse_ce(y ~ s1 + s2 + s3 + s4 + s5 - 1, data = dat)

  # Rebuild the dual gradient g(lambda) = X p(lambda) - y and finite-difference
  logq <- log(rep(1 / 5, 5))
  probs <- function(l) { e <- exp(logq + as.vector(l %*% X)); e / sum(e) }
  g     <- function(l) as.vector(X %*% probs(l) - dat$y)
  lam <- fit$lambda_hat; h <- 1e-6
  Hfd <- sapply(seq_along(lam), function(j) {
    lp <- lam; lp[j] <- lp[j] + h; (g(lp) - g(lam)) / h
  })
  Hfd <- (Hfd + t(Hfd)) / 2
  expect_equal(unname(fit$hessian), unname(Hfd), tolerance = 1e-3)
})

test_that("vcov equals the inverse of the information matrix (full rank)", {
  X <- rbind(c(1, 2, 3, 4, 5),
             c(1, 4, 9, 16, 25))
  p_true <- exp(as.vector(crossprod(X, c(0.2, -0.05))))
  p_true <- p_true / sum(p_true)
  dat <- data.frame(y  = as.vector(X %*% p_true),
                    s1 = X[, 1], s2 = X[, 2], s3 = X[, 3],
                    s4 = X[, 4], s5 = X[, 5])
  fit <- inverse_ce(y ~ s1 + s2 + s3 + s4 + s5 - 1, data = dat)

  V <- vcov(fit)
  expect_equal(dim(V), c(2L, 2L))
  expect_equal(unname(V), unname(solve(fit$hessian)), tolerance = 1e-8)
  expect_equal(V, t(V), tolerance = 1e-10)
})

test_that("information-matrix SEs of lambda and p are finite and positive", {
  X <- rbind(c(1, 2, 3, 4, 5),
             c(1, 4, 9, 16, 25))
  p_true <- exp(as.vector(crossprod(X, c(0.2, -0.05))))
  p_true <- p_true / sum(p_true)
  dat <- data.frame(y  = as.vector(X %*% p_true),
                    s1 = X[, 1], s2 = X[, 2], s3 = X[, 3],
                    s4 = X[, 4], s5 = X[, 5])
  fit <- inverse_ce(y ~ s1 + s2 + s3 + s4 + s5 - 1, data = dat)

  expect_length(fit$se_lambda, 2)
  expect_length(fit$se_p, 5)
  expect_true(all(is.finite(fit$se_lambda) & fit$se_lambda > 0))
  expect_true(all(is.finite(fit$se_p) & fit$se_p >= 0))
  # se_lambda = sqrt(diag(vcov_lambda))
  expect_equal(unname(fit$se_lambda),
               unname(sqrt(diag(fit$vcov_lambda))), tolerance = 1e-10)
})

test_that("a more informative moment (larger Var_p(x_t)) has a smaller SE(lambda_t)", {
  # Row 2 (squares) has far larger p-weighted variance than row 1 (linear),
  # so lambda_2 is more sharply identified: SE(lambda_2) < SE(lambda_1).
  X <- rbind(c(1, 2, 3, 4, 5),
             c(1, 4, 9, 16, 25))
  p_true <- exp(as.vector(crossprod(X, c(0.2, -0.05))))
  p_true <- p_true / sum(p_true)
  dat <- data.frame(y  = as.vector(X %*% p_true),
                    s1 = X[, 1], s2 = X[, 2], s3 = X[, 3],
                    s4 = X[, 4], s5 = X[, 5])
  fit <- inverse_ce(y ~ s1 + s2 + s3 + s4 + s5 - 1, data = dat)

  Ep   <- as.vector(X %*% fit$p_hat)
  varx <- vapply(1:2, function(t) sum(fit$p_hat * (X[t, ] - Ep[t])^2), numeric(1))
  expect_equal(which.max(varx), unname(which.min(fit$se_lambda)))
})

test_that("delta-method se_p uses the correct Jacobian dp/dlambda", {
  X <- rbind(c(1, 2, 3, 4, 5),
             c(1, 4, 9, 16, 25))
  p_true <- exp(as.vector(crossprod(X, c(0.2, -0.05))))
  p_true <- p_true / sum(p_true)
  dat <- data.frame(y  = as.vector(X %*% p_true),
                    s1 = X[, 1], s2 = X[, 2], s3 = X[, 3],
                    s4 = X[, 4], s5 = X[, 5])
  fit <- inverse_ce(y ~ s1 + s2 + s3 + s4 + s5 - 1, data = dat)

  # analytic J_{k,t} = p_k (x_tk - E_p[x_t]) vs finite difference of p(lambda)
  logq <- log(rep(1 / 5, 5))
  probs <- function(l) { e <- exp(logq + as.vector(l %*% X)); e / sum(e) }
  Ep <- as.vector(X %*% fit$p_hat)
  J  <- fit$p_hat * (t(X) - matrix(Ep, 5, 2, byrow = TRUE))
  lam <- fit$lambda_hat; h <- 1e-6
  Jfd <- sapply(1:2, function(j) { lp <- lam; lp[j] <- lp[j] + h;
                                   (probs(lp) - probs(lam)) / h })
  expect_equal(unname(J), unname(Jfd), tolerance = 1e-4)
  # and se_p = sqrt(diag(J vcov J'))
  se_p_check <- sqrt(pmax(diag(J %*% fit$vcov_lambda %*% t(J)), 0))
  expect_equal(unname(fit$se_p), unname(se_p_check), tolerance = 1e-10)
})

test_that("singular information matrix (T >= K) yields NA SEs and NULL vcov, no error", {
  # T = 3 moments, K = 3 states: Cov_p(moments) has rank <= K - 1 = 2 < 3.
  X <- rbind(c(1, 2, 3),
             c(1, 4, 9),
             c(1, 8, 27))
  dat <- data.frame(y = c(2, 4, 8), a = X[, 1], b = X[, 2], c = X[, 3])
  fit <- suppressWarnings(inverse_ce(y ~ a + b + c - 1, data = dat))

  expect_true(all(is.na(fit$se_lambda)))
  expect_true(all(is.na(fit$se_p)))
  expect_null(fit$vcov_lambda)
  expect_null(vcov(fit))
})

test_that("a constant moment row also makes the information matrix singular", {
  # Moment row 1 is constant across states -> zero p-variance -> singular I.
  X <- rbind(c(1, 1, 1, 1, 1),
             c(1, 2, 3, 4, 5))
  dat <- data.frame(y = c(1, 3), a = X[, 1], b = X[, 2], c = X[, 3],
                    d = X[, 4], e = X[, 5])
  fit <- suppressWarnings(inverse_ce(y ~ a + b + c + d + e - 1, data = dat))
  expect_true(all(is.na(fit$se_lambda)))
  expect_null(vcov(fit))
})

test_that("affinely-dependent moment rows make the information matrix singular", {
  # Row 2 = 6 - Row 1, so the two moments are affinely dependent: Cov_p is
  # rank 1 (not 2) and I(lambda) is singular -> NA SEs, NULL vcov.
  X <- rbind(c(1, 2, 3, 4, 5),
             c(5, 4, 3, 2, 1))
  p_true <- c(.30, .25, .20, .15, .10)
  dat <- data.frame(y = as.vector(X %*% p_true),
                    x1 = X[, 1], x2 = X[, 2], x3 = X[, 3],
                    x4 = X[, 4], x5 = X[, 5])
  fit <- inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat)
  expect_true(all(is.na(fit$se_lambda)))
  expect_null(vcov(fit))
})

test_that("normalized entropy S lies in [0, 1]", {
  X <- matrix(c(1, 2, 3, 4, 5,
                5, 4, 3, 2, 1), nrow = 2, byrow = TRUE)
  p_true <- c(.30, .25, .20, .15, .10)
  dat <- data.frame(y = as.vector(X %*% p_true),
                    x1 = X[, 1], x2 = X[, 2], x3 = X[, 3],
                    x4 = X[, 4], x5 = X[, 5])
  fit <- inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat)

  expect_gte(fit$S, 0)
  expect_lte(fit$S, 1)
})

test_that("fano_bounds returns a one-row frame with a valid Fano bound", {
  X <- rbind(c(1, 2, 3, 4, 5),
             c(1, 4, 9, 16, 25))
  p_true <- exp(as.vector(crossprod(X, c(0.2, -0.05))))
  p_true <- p_true / sum(p_true)
  dat <- data.frame(y  = as.vector(X %*% p_true),
                    s1 = X[, 1], s2 = X[, 2], s3 = X[, 3],
                    s4 = X[, 4], s5 = X[, 5])
  fit <- inverse_ce(y ~ s1 + s2 + s3 + s4 + s5 - 1, data = dat)

  fb <- fano_bounds(fit)
  expect_s3_class(fb, "data.frame")
  expect_equal(nrow(fb), 1L)
  expect_true(all(c("p_max", "pe", "H", "S", "pe_lower") %in% names(fb)))

  # weak Fano: pe >= S(p) - log2/logK  (the reported pe_lower, clamped at 0)
  K <- length(coef(fit))
  expect_gte(fb$pe, fb$pe_lower - 1e-9)
  expect_equal(fb$pe_lower, max(0, fb$S - log(2) / log(K)), tolerance = 1e-10)
  expect_equal(fb$pe, 1 - max(coef(fit)), tolerance = 1e-10)

  # strong Fano: H(pe) + pe log(K-1) >= H(p)
  pe  <- fb$pe
  Hpe <- if (pe > 0 && pe < 1) -pe * log(pe) - (1 - pe) * log(1 - pe) else 0
  expect_gte(Hpe + pe * log(K - 1), fb$H - 1e-9)

  # "overall" attribute mirrors the single row
  ov <- attr(fb, "overall")
  expect_equal(unname(ov["mean_pe"]), fb$pe, tolerance = 1e-12)
  expect_equal(unname(ov["mean_pe_lower"]), fb$pe_lower, tolerance = 1e-12)
})

test_that("input validation throws informative errors", {
  X <- matrix(c(1, 2, 3, 4, 5,
                5, 4, 3, 2, 1), nrow = 2, byrow = TRUE)
  dat <- data.frame(y = as.vector(X %*% rep(1 / 5, 5)),
                    x1 = X[, 1], x2 = X[, 2], x3 = X[, 3],
                    x4 = X[, 4], x5 = X[, 5])

  # wrong-length prior
  expect_error(
    inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat, p0 = rep(1 / 3, 3)),
    "must equal the number of states"
  )
  # non-positive prior mass
  expect_error(
    inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat,
               p0 = c(.4, .3, .2, .1, 0)),
    "strictly positive"
  )
})

test_that("intercept in formula triggers a warning", {
  # Single moment, K = 4 states; with the auto-intercept K = 5 > T + 1 = 2,
  # so the only warning is the intended intercept message.
  dat <- data.frame(y = 5, x1 = 2, x2 = 4, x3 = 6, x4 = 8)
  expect_warning(
    inverse_ce(y ~ x1 + x2 + x3 + x4, data = dat),
    "intercept"
  )
})

test_that("non-uniform prior that does not sum to 1 is renormalized", {
  X <- matrix(c(1, 2, 3, 4, 5), nrow = 1)
  dat <- data.frame(y = 3, x1 = 1, x2 = 2, x3 = 3, x4 = 4, x5 = 5)
  expect_warning(
    fit <- inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat,
                      p0 = c(2, 2, 2, 2, 2)),
    "renormaliz"
  )
  expect_equal(sum(fit$p0), 1, tolerance = 1e-12)
})

test_that("S3 methods run and return correctly shaped objects", {
  # Non-collinear moment rows so I(lambda) is full rank (vcov is 2 x 2).
  X <- rbind(c(1, 2, 3, 4, 5),
             c(1, 4, 9, 16, 25))
  p_true <- exp(as.vector(crossprod(X, c(0.2, -0.05))))
  p_true <- p_true / sum(p_true)
  dat <- data.frame(y = as.vector(X %*% p_true),
                    x1 = X[, 1], x2 = X[, 2], x3 = X[, 3],
                    x4 = X[, 4], x5 = X[, 5])
  fit <- inverse_ce(y ~ x1 + x2 + x3 + x4 + x5 - 1, data = dat)

  expect_length(coef(fit), 5)
  expect_length(fitted(fit), 2)
  expect_length(residuals(fit), 2)
  expect_equal(dim(vcov(fit)), c(2L, 2L))
  expect_output(print(fit), "Estimated probabilities")
  s <- summary(fit)
  expect_s3_class(s, "summary.inverse_ce")
  expect_output(print(s), "Coefficients")
  expect_output(print(s), "Std. Error")     # SE columns printed
  expect_output(print(s), "Fano")           # Fano line printed
})
