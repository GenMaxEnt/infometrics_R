# Tests for markov_gce() — GCE (noisy-moment) Markov transition matrix from a
# panel with covariates (Golan 2008, Sec. 7.7.1, eqs. 7.25a / 7.28).

# Vectorized balanced share-panel with time-varying covariates.
sim_gce_panel <- function(P, N, Tn, seed) {
  set.seed(seed)
  K <- nrow(P)
  Y <- matrix(rgamma(N * K, 1), N, K); Y <- Y / rowSums(Y); inc <- rnorm(N)
  id <- rep(seq_len(N), each = Tn); tm <- rep(seq_len(Tn), N)
  sh <- matrix(0, N * Tn, K); income <- numeric(N * Tn); infl <- numeric(N * Tn)
  for (t in seq_len(Tn)) {
    idx <- which(tm == t); inc <- 0.8 * inc + rnorm(N, 0, 0.5)
    sh[idx, ] <- Y; income[idx] <- inc; infl[idx] <- 0.02 * t
    Y <- (Y %*% P) * exp(matrix(rnorm(N * K, 0, 0.03), N, K)); Y <- Y / rowSums(Y)
  }
  d <- data.frame(i = id, t = tm)
  for (k in seq_len(K)) d[[paste0("y", k)]] <- sh[, k]
  d$income <- income; d$infl <- infl
  d
}

P_true <- matrix(c(.6, .3, .1, .2, .5, .3, .1, .3, .6), 3, byrow = TRUE)
st3 <- c("y1", "y2", "y3")
cv  <- c("income", "infl")

test_that("converges and satisfies the noisy first-order conditions", {
  panel <- sim_gce_panel(P_true, N = 120L, Tn = 5L, seed = 1)
  fit <- markov_gce(panel, "i", "t", st3, covariates = cv,
                    v = c(-0.2, 0, 0.2), nu = 0.2)

  expect_s3_class(fit, "markov_gce")
  expect_s3_class(fit, "infometrics")
  expect_equal(fit$method, "dual")
  expect_true(fit$converged)
  expect_lt(max(abs(rowSums(coef(fit)) - 1)), 1e-8)     # row-stochastic
  expect_lt(fit$foc_residual, 1e-3)                     # noisy FOC satisfied
  expect_gt(fit$moment_residual, 0)                     # exact-moment gap expected nonzero
  expect_equal(dim(fit$epsilon), c(fit$n_transitions, 3L))
  expect_equal(dim(fit$lambda), c(2L, 3L))              # S x K
})

test_that("larger nu pulls the signal toward the prior (higher S)", {
  panel <- sim_gce_panel(P_true, N = 120L, Tn = 5L, seed = 2)
  f_lo <- markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-0.2, 0, 0.2), nu = 0.2)
  f_hi <- markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-0.2, 0, 0.2), nu = 0.8)
  expect_gt(f_hi$S, f_lo$S)
})

test_that("a wider noise support pulls the signal toward the prior (higher S)", {
  panel <- sim_gce_panel(P_true, N = 120L, Tn = 5L, seed = 3)
  f_narrow <- markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-0.2, 0, 0.2), nu = 0.5)
  f_wide   <- markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-1, 0, 1),     nu = 0.5)
  expect_gt(f_wide$S, f_narrow$S)
})

test_that("canonical aliases are present and consistent", {
  panel <- sim_gce_panel(P_true, N = 100L, Tn = 5L, seed = 4)
  fit <- markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-0.2, 0, 0.2), nu = 0.3)
  expect_identical(fit$p_hat, fit$p)
  expect_identical(fit$lambda_hat, fit$lambda)
  expect_equal(fit$H_signal, fit$entropy)
  expect_equal(fit$H_signal, fit$H_p)
})

test_that("input validation throws informative errors", {
  panel <- sim_gce_panel(P_true, N = 80L, Tn = 4L, seed = 6)

  expect_error(markov_gce(panel, "i", "t", st3), "covariates' is required")
  expect_error(markov_gce(panel, "i", "t", st3, covariates = cv, nu = 0),
               "in \\(0, 1\\)")
  expect_error(markov_gce(panel, "i", "t", st3, covariates = cv, nu = 1),
               "in \\(0, 1\\)")
  expect_error(markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-1, 0, 2)),
               "symmetric")
  expect_error(markov_gce(panel, "i", "t", st3, covariates = cv, v = 0.5),
               "at least two support points")
  expect_error(markov_gce(panel, "i", "t", st3, covariates = cv,
                          p0 = matrix(1 / 2, 2, 2)), "p0 must be K x K")
  expect_error(markov_gce(panel, "i", "t", st3, covariates = cv,
                          v = c(-1, 0, 1), w0 = c(.5, .5)), "length\\(v\\)")
})

test_that("S3 methods run and return correctly shaped objects", {
  panel <- sim_gce_panel(P_true, N = 100L, Tn = 5L, seed = 7)
  fit <- markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-0.2, 0, 0.2), nu = 0.3)
  nt <- fit$n_transitions

  expect_equal(dim(coef(fit)), c(3L, 3L))               # K x K
  expect_equal(dim(fitted(fit)), c(nt, 3L))
  expect_equal(dim(residuals(fit)), c(nt, 3L))
  expect_output(print(fit), "Generalized Cross-Entropy")
  s <- summary(fit)
  expect_s3_class(s, "summary.markov_gce")
  expect_output(print(s), "GCE transition matrix")
})

# ---- fano_bounds() and margins() with sandwich SEs --------------------------

Hbin <- function(pe) ifelse(pe > 0 & pe < 1, -pe * log(pe) - (1 - pe) * log(1 - pe), 0)

test_that("fano_bounds returns per-origin bounds that hold", {
  panel <- sim_gce_panel(P_true, N = 120L, Tn = 5L, seed = 1)
  fit <- markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-0.2, 0, 0.2), nu = 0.2)
  fb <- fano_bounds(fit)
  expect_equal(nrow(fb), length(st3))
  expect_true(all(fb$pe >= fb$pe_lower - 1e-9))
  expect_true(all(Hbin(fb$pe) + fb$pe * log(length(st3) - 1) >= fb$H - 1e-9))
  expect_named(attr(fb, "overall"), c("mean_pe", "mean_pe_lower", "S_system"))
})

test_that("margins point estimates match the softmax-Jacobian identity", {
  panel <- sim_gce_panel(P_true, N = 150L, Tn = 6L, seed = 2)
  fit <- markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-0.2, 0, 0.2), nu = 0.3)
  ME <- margins(fit, average = FALSE)
  expect_equal(dim(ME), c(3L, 3L, length(cv)))
  P <- fit$p_hat; lam <- fit$lambda; nu <- fit$nu
  for (k in 1:3) { Jk <- diag(P[k, ]) - outer(P[k, ], P[k, ])
    for (s in seq_along(cv)) expect_equal(unname(ME[k, , s]),
                                          unname(as.vector(Jk %*% lam[s, ]) / nu), tolerance = 1e-10) }
  expect_equal(dim(margins(fit)), c(length(cv), 3L))            # AME S x K
})

test_that("margins sandwich SEs: present, finite, and structured", {
  panel <- sim_gce_panel(P_true, N = 150L, Tn = 6L, seed = 2)
  fit <- markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-0.2, 0, 0.2), nu = 0.3)
  m <- margins(fit, se = TRUE)
  expect_s3_class(m, "margins_gce")
  expect_equal(m$se_method, "sandwich")
  expect_equal(m$estimate, margins(fit), tolerance = 1e-10)
  expect_equal(dim(m$se), c(length(cv), 3L))
  expect_true(all(is.finite(m$se) & m$se >= 0))
  expect_equal(m$z, m$estimate / m$se, tolerance = 1e-12)
  expect_output(print(m), "marginal effects")
  expect_error(margins(fit, se = TRUE, average = FALSE), "average = TRUE")
})

# well-identified 3-covariate panel (all AR(1)) for the SE-comparison tests;
# avoids the deterministic 'infl' column, whose weak identification makes the
# sandwich under-state its SE (a documented limitation, not a bug).
sim_gce_panel3 <- function(P, N, Tn, seed) {
  set.seed(seed); K <- nrow(P)
  Y <- matrix(rgamma(N * K, 1), N, K); Y <- Y / rowSums(Y)
  id <- rep(seq_len(N), each = Tn); tm <- rep(seq_len(Tn), N)
  sh <- matrix(0, N * Tn, K); Z <- matrix(0, N * Tn, 3); C <- matrix(rnorm(N * 3), N, 3)
  for (t in seq_len(Tn)) {
    idx <- which(tm == t); C <- 0.7 * C + matrix(rnorm(N * 3, 0, 0.6), N, 3)
    sh[idx, ] <- Y; Z[idx, ] <- C
    Y <- (Y %*% P) * exp(matrix(rnorm(N * K, 0, 0.03), N, K)); Y <- Y / rowSums(Y)
  }
  d <- data.frame(i = id, t = tm)
  for (k in seq_len(K)) d[[paste0("y", k)]] <- sh[, k]
  for (s in 1:3) d[[paste0("z", s)]] <- Z[, s]
  d
}
cv3 <- c("z1", "z2", "z3")

test_that("bootstrap se_method runs and matches the sandwich in aggregate", {
  panel <- sim_gce_panel3(P_true, N = 150L, Tn = 6L, seed = 2)
  fit <- markov_gce(panel, "i", "t", st3, covariates = cv3, v = c(-0.2, 0, 0.2), nu = 0.3)
  ms <- margins(fit, se = TRUE)                         # sandwich (default)
  set.seed(7)
  mb <- margins(fit, se = TRUE, se_method = "bootstrap", B = 80L)
  expect_s3_class(mb, "margins_gce")
  expect_equal(mb$se_method, "bootstrap")
  expect_equal(dim(mb$se), c(length(cv3), 3L))
  expect_true(all(is.finite(mb$se) & mb$se >= 0))
  expect_true(mb$n_boot >= 2L)
  expect_equal(mb$estimate, ms$estimate, tolerance = 1e-10)   # same point estimate
  # the internal unit-block bootstrap agrees with the sandwich in aggregate
  expect_true(median(ms$se / mb$se) > 0.5 && median(ms$se / mb$se) < 2)
  expect_output(print(mb), "bootstrap resamples used")
})

test_that("bootstrap se_method agrees with a public refit-from-data bootstrap", {
  panel <- sim_gce_panel3(P_true, N = 150L, Tn = 6L, seed = 2)
  fit <- markov_gce(panel, "i", "t", st3, covariates = cv3, v = c(-0.2, 0, 0.2), nu = 0.3)
  set.seed(7)
  se_internal <- margins(fit, se = TRUE, se_method = "bootstrap", B = 80L)$se
  # reference: resample the raw panel by unit and refit the public estimator
  units <- unique(panel$i); B <- 80L; AMEb <- array(NA_real_, c(length(cv3), 3L, B))
  set.seed(7)
  for (b in seq_len(B)) {
    bu <- sample(units, length(units), replace = TRUE)
    dl <- do.call(rbind, lapply(seq_along(bu), function(m) { d <- panel[panel$i == bu[m], ]; d$i <- m; d }))
    fb <- tryCatch(markov_gce(dl, "i", "t", st3, covariates = cv3, v = c(-0.2, 0, 0.2), nu = 0.3),
                   error = function(e) NULL)
    if (!is.null(fb) && fb$converged) AMEb[, , b] <- margins(fb)
  }
  se_public <- apply(AMEb, c(1, 2), function(z) sd(z, na.rm = TRUE))
  expect_true(median(se_internal / se_public) > 0.6 && median(se_internal / se_public) < 1.6)
})

# ---- analytic-Hessian standard errors of lambda (Golan 2008 eqs. 4.6/4.7) ----

test_that("analytic Hessian (signal + noise) equals the finite-difference Hessian", {
  panel <- sim_gce_panel(P_true, N = 150L, Tn = 6L, seed = 1)
  fit <- markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-0.2, 0, 0.2), nu = 0.3)
  sd <- fit$se_data; lhat <- as.vector(fit$lambda * sd$cov_scale)
  Han <- infometrics:::.markov_info(sd, lhat)
  grad <- function(l) { pc <- infometrics:::.markov_P_eps(l, sd)
    as.vector(sd$M - crossprod(sd$A, pc$P) - crossprod(sd$Z1s, pc$eps)) }
  Hfd <- -infometrics:::.fd_jac(grad, lhat, h = 1e-6)      # tighter step: finite-diff error ~1e-7
  expect_lt(max(abs(Han - Hfd)) / max(abs(Hfd)), 1e-5)     # analytic exact up to finite-diff error
})

test_that("se_lambda is full S x K (no reference NA), finite, and PD vcov", {
  panel <- sim_gce_panel(P_true, N = 150L, Tn = 6L, seed = 1)
  fit <- markov_gce(panel, "i", "t", st3, covariates = cv, v = c(-0.2, 0, 0.2), nu = 0.3)
  expect_equal(dim(fit$se_lambda), c(length(cv), 3L))
  expect_true(all(is.finite(fit$se_lambda) & fit$se_lambda >= 0))   # no NA (identified)
  expect_true(all(eigen(fit$vcov_lambda, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_output(print(summary(fit)), "analytic-Hessian unit-clustered")
})

test_that("clustered analytic se_lambda matches a unit bootstrap; non-clustered overstates", {
  panel <- sim_gce_panel3(P_true, N = 150L, Tn = 6L, seed = 2)   # well-identified covariates
  fit <- markov_gce(panel, "i", "t", st3, covariates = cv3, v = c(-0.2, 0, 0.2), nu = 0.3)
  units <- unique(panel$i); B <- 120L; LB <- array(NA_real_, c(length(cv3), 3L, B))
  set.seed(4)
  for (b in seq_len(B)) {
    bu <- sample(units, length(units), replace = TRUE)
    dl <- do.call(rbind, lapply(seq_along(bu), function(m) { d <- panel[panel$i == bu[m], ]; d$i <- m; d }))
    fb <- tryCatch(markov_gce(dl, "i", "t", st3, covariates = cv3, v = c(-0.2, 0, 0.2), nu = 0.3),
                   error = function(e) NULL)
    if (!is.null(fb) && fb$converged) LB[, , b] <- fb$lambda
  }
  seb <- apply(LB, c(1, 2), function(z) sd(z, na.rm = TRUE))
  expect_true(median(fit$se_lambda / seb) > 0.7 && median(fit$se_lambda / seb) < 1.4)  # clustered ~ boot

  # a NON-clustered analytic sandwich overstates (documents why clustering is used)
  sd <- fit$se_data; S <- length(cv3); Kk <- 3L; lhat <- as.vector(fit$lambda * sd$cov_scale)
  H <- infometrics:::.markov_info(sd, lhat); Hi <- solve(H)
  pc <- infometrics:::.markov_P_eps(lhat, sd); fitm <- sd$XX %*% pc$P
  G <- matrix(0, nrow(sd$XX), S * Kk)
  for (j in 1:Kk) G[, ((j - 1) * S + 1):(j * S)] <- sd$Z1s * (sd$YY[, j] - pc$eps[, j]) - sd$Z2s * fitm[, j]
  D <- 1 / sd$cov_scale[rep(1:S, Kk)]
  se_nc <- matrix(sqrt(diag(outer(D, D) * (Hi %*% crossprod(G) %*% Hi))), S, Kk)
  expect_gt(median(se_nc / fit$se_lambda), 1.1)            # non-clustered > clustered
})
