# Tests for mixed_gce() — doubly-reparameterized GCE mixed model.

# Self-consistent DGP: pick lambda*, root-find rho* so sum_j p_ij = 1, build the
# implied signal shares Pstar, and set Y = Pstar (noise-free). A well-specified
# fit then recovers lambda*/Pstar in the shrink-to-pure (u -> 0) limit.
sim_mixed <- function(N = 40L, J = 3L, K = 2L, nu = 0.5,
                      s = c(0, 0.5, 1), seed = 1) {
  set.seed(seed)
  X <- array(0, dim = c(N, J, K))
  X[, , 1] <- 1
  for (k in 2:K) X[, , k] <- matrix(rnorm(N * J), N, J)
  lam <- matrix(seq(-0.4, 0.5, length.out = K * J), K, J)
  V <- matrix(0, N, J)
  for (j in seq_len(J)) V[, j] <- X[, j, ] %*% lam[, j]
  pfun <- function(r, Vi) sapply(seq_len(J), function(j) {
    a <- s * (Vi[j] + r) / nu; a <- a - max(a); th <- exp(a) / sum(exp(a)); sum(s * th)
  })
  Pstar <- matrix(0, N, J)
  for (i in seq_len(N)) {
    ri <- stats::uniroot(function(r) sum(pfun(r, V[i, ])) - 1, c(-60, 60))$root
    Pstar[i, ] <- pfun(ri, V[i, ])
  }
  list(Y = Pstar, X = X, lam = lam, Pstar = Pstar, N = N, J = J, K = K, s = s, nu = nu)
}

moment_resid <- function(fit) {
  R <- fit$y_mat - fit$p_hat - fit$e_hat
  m <- numeric(0)
  for (j in seq_len(fit$J)) m <- c(m, colSums(matrix(fit$X[, j, ], nrow = fit$N) * R[, j]))
  max(abs(m))
}

test_that("converges and satisfies the dual first-order conditions", {
  d   <- sim_mixed()
  fit <- mixed_gce(d$Y, d$X, nu = d$nu)
  expect_s3_class(fit, "mixed_gce")
  expect_s3_class(fit, "infometrics")
  expect_equal(fit$method, "dual")
  expect_true(fit$converged)
  g <- infometrics:::.mixed_gce_grad(c(fit$rho, as.vector(fit$lambda)),
                                     d$Y, d$X, fit$N, fit$J, fit$K, fit$s, fit$u,
                                     fit$theta0, fit$w0, fit$nu)
  expect_lt(max(abs(g)), 1e-4)
})

test_that("normalization, probability validity, and distribution fibers hold", {
  d   <- sim_mixed()
  fit <- mixed_gce(d$Y, d$X, nu = d$nu)
  expect_lt(max(abs(rowSums(fit$p_hat) - 1)), 1e-6)          # sum_j p_ij = 1
  expect_true(all(fit$p_hat >= -1e-9 & fit$p_hat <= 1 + 1e-9))
  expect_lt(max(abs(apply(fit$theta_hat, c(1, 2), sum) - 1)), 1e-8)
  expect_lt(max(abs(apply(fit$w_hat, c(1, 2), sum) - 1)), 1e-8)
})

test_that("satisfies the data moment condition", {
  d   <- sim_mixed()
  fit <- mixed_gce(d$Y, d$X, nu = d$nu)
  expect_lt(moment_resid(fit), 1e-4)
})

test_that("recovers lambda*/Pstar in the shrink-to-pure limit", {
  d  <- sim_mixed()
  f1 <- mixed_gce(d$Y, d$X, u = c(-0.10, 0, 0.10), nu = d$nu)
  f2 <- mixed_gce(d$Y, d$X, u = c(-0.01, 0, 0.01), nu = d$nu)
  # tighter noise support -> closer to the pure-signal DGP
  expect_lt(max(abs(f2$p_hat - d$Pstar)), 1e-2)
  expect_lt(max(abs(f2$lambda[2, ] - d$lam[2, ])), 1e-2)      # covariate row
  expect_lt(max(abs(f2$p_hat - d$Pstar)), max(abs(f1$p_hat - d$Pstar)))
})

test_that("GME equals GCE with uniform priors; a non-uniform prior shifts p", {
  d  <- sim_mixed()
  u0 <- c(-0.05, 0, 0.05)
  f_gme <- mixed_gce(d$Y, d$X, u = u0, nu = d$nu)
  th0   <- array(0, dim = c(d$N, d$J, 3)); for (m in 1:3) th0[, , m] <- c(0.6, 0.3, 0.1)[m]
  th0u  <- array(1 / 3, dim = c(d$N, d$J, 3))
  f_uni <- mixed_gce(d$Y, d$X, theta0 = th0u, u = u0, nu = d$nu)
  f_gce <- mixed_gce(d$Y, d$X, theta0 = th0,  u = u0, nu = d$nu)
  expect_equal(f_gme$p_hat, f_uni$p_hat, tolerance = 1e-6)    # uniform prior == GME
  expect_gt(max(abs(f_gce$p_hat - f_gme$p_hat)), 1e-6)        # informative prior shifts p
})

test_that("flexible s/u supports: NULL, count, and explicit vector", {
  d   <- sim_mixed()
  f_default <- mixed_gce(d$Y, d$X, nu = d$nu)
  f_count   <- mixed_gce(d$Y, d$X, u = 5L, s = 4L, nu = d$nu)
  f_vec     <- mixed_gce(d$Y, d$X, u = c(-0.5, 0, 0.5), nu = d$nu)
  expect_equal(f_default$H, 3L); expect_equal(f_default$M, 3L)
  expect_equal(f_count$H, 5L);   expect_equal(f_count$M, 4L)
  expect_equal(f_vec$u, c(-0.5, 0, 0.5))
})

test_that("standard errors are present, PD, and finite", {
  d   <- sim_mixed()
  fit <- mixed_gce(d$Y, d$X, nu = d$nu)
  V <- vcov(fit)
  expect_equal(dim(V), rep(fit$N + fit$K * fit$J, 2L))
  expect_equal(V, t(V))
  expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_equal(dim(fit$se_lambda), c(fit$K, fit$J))
  expect_true(all(is.finite(fit$se_lambda) & fit$se_lambda >= 0))
  expect_length(fit$se_rho, fit$N)
  expect_true(all(is.finite(fit$se_rho) & fit$se_rho >= 0))
})

test_that("margins reproduce a central-difference derivative", {
  d   <- sim_mixed(N = 25L)
  fit <- mixed_gce(d$Y, d$X, nu = d$nu)
  ME  <- margins(fit, average = FALSE)
  expect_equal(dim(ME), c(fit$N, fit$J, fit$K))
  expect_equal(dim(margins(fit)), c(fit$K, fit$J))
  lam <- fit$lambda; rho <- fit$rho
  p_of_x <- function(i, j, k, delta) {
    Xp <- d$X; Xp[i, j, k] <- Xp[i, j, k] + delta
    Vp <- matrix(0, fit$N, fit$J)
    for (jj in seq_len(fit$J)) Vp[, jj] <- Xp[, jj, ] %*% lam[, jj]
    infometrics:::.mixed_gce_signal(Vp + rho, fit$s, fit$theta0, fit$nu)$p[i, j]
  }
  h <- 1e-6; err <- 0
  for (t in 1:6) {
    i <- ((t * 7L) %% fit$N) + 1L; j <- (t %% fit$J) + 1L; k <- (t %% fit$K) + 1L
    num <- (p_of_x(i, j, k, h) - p_of_x(i, j, k, -h)) / (2 * h)
    err <- max(err, abs(num - ME[i, j, k]))
  }
  expect_lt(err, 1e-5)
})

test_that("input validation throws informative errors", {
  d <- sim_mixed(N = 20L)
  expect_error(mixed_gce(d$Y[-1, ], d$X, nu = 0.5), "dim\\(Y\\)")
  expect_error(mixed_gce(d$Y, d$X[, , 1], nu = 0.5), "3-dimensional")
  expect_error(mixed_gce(d$Y, d$X, nu = 0), "inside \\(0, 1\\)")
  expect_error(mixed_gce(d$Y, d$X, nu = 1), "inside \\(0, 1\\)")
  expect_error(mixed_gce(d$Y, d$X, s = c(0, 0.5, 1.5)), "\\[0, 1\\]")
  expect_error(mixed_gce(d$Y, d$X, u = c(-1, 0, 0.5)), "symmetric")
  expect_error(mixed_gce(d$Y, d$X, theta0 = array(1 / 2, dim = c(d$N, d$J, 2))),
               "dim\\(theta0\\)")
})

test_that("canonical aliases and S3 methods are consistent", {
  d   <- sim_mixed(N = 20L)
  fit <- mixed_gce(d$Y, d$X, nu = d$nu)
  expect_identical(coef(fit), fit$lambda)
  expect_identical(fit$lambda_hat, fit$lambda)
  expect_identical(fit$p_hat, fit$p)
  expect_identical(fit$objective, fit$value)
  expect_equal(fit$converged, fit$convergence == 0L)
  expect_length(fit$H_signal, fit$N)
  expect_equal(dim(fitted(fit)), c(fit$N, fit$J))
  expect_equal(dim(residuals(fit)), c(fit$N, fit$J))
  expect_equal(unname(fitted(fit) + residuals(fit) + fit$e_hat), unname(d$Y),
               tolerance = 1e-10)
  expect_output(print(fit), "Mixed")
  expect_output(print(summary(fit)), "Category")
})

# ---- Fano error bounds (Golan 2008, sec 7.5) --------------------------------

test_that("fano_bounds returns per-observation bounds that hold", {
  d   <- sim_mixed()
  fit <- mixed_gce(d$Y, d$X, nu = d$nu)
  fb  <- fano_bounds(fit)

  expect_s3_class(fb, "data.frame")
  expect_equal(nrow(fb), fit$N)
  expect_named(fb, c("p_max", "pe", "H", "S", "pe_lower"))
  expect_equal(fb$pe, 1 - apply(fit$p_hat, 1, max), tolerance = 1e-12)
  expect_true(all(fb$S >= 0 & fb$S <= 1 + 1e-12))
  expect_true(all(fb$pe >= fb$pe_lower - 1e-9))                 # weak Fano bound
  Hbin <- function(pe) ifelse(pe > 0 & pe < 1,
                              -pe * log(pe) - (1 - pe) * log(1 - pe), 0)
  expect_true(all(Hbin(fb$pe) + fb$pe * log(fit$J - 1) >= fb$H - 1e-9))  # strong
  expect_named(attr(fb, "overall"), c("mean_pe", "mean_pe_lower", "S_system"))
})

# ---- marginal-effect standard errors (robust sandwich) ----------------------

test_that("margins(se = TRUE) sandwich: estimates, SEs, and structure", {
  d   <- sim_mixed()
  fit <- mixed_gce(d$Y, d$X, nu = d$nu)
  m   <- margins(fit, se = TRUE)                       # default se_method = "sandwich"

  expect_s3_class(m, "margins_gce")
  expect_equal(m$se_method, "sandwich")
  expect_equal(m$estimate, margins(fit), tolerance = 1e-10)
  expect_equal(dim(m$se), c(fit$K, fit$J))
  expect_true(all(is.finite(m$se) & m$se >= 0))
  expect_equal(m$z, m$estimate / m$se, tolerance = 1e-12)
  expect_output(print(m), "marginal effects")

  # the sandwich SE is far below the naive Hessian delta (which is ~7-9x too big)
  V_naive <- fit$vcov                                  # full Hessian-inverse
  par_hat <- c(fit$rho, as.vector(fit$lambda)); np <- length(par_hat); hh <- 1e-5
  Jac <- matrix(0, fit$K * fit$J, np)
  for (mm in seq_len(np)) { e <- numeric(np); e[mm] <- hh
    Jac[, mm] <- (infometrics:::.mixed_gce_me_vec(par_hat + e, fit$X, fit$s, fit$theta0, fit$nu, fit$N, fit$J, fit$K) -
                  infometrics:::.mixed_gce_me_vec(par_hat - e, fit$X, fit$s, fit$theta0, fit$nu, fit$N, fit$J, fit$K)) / (2 * hh)
  }
  se_naive <- matrix(sqrt(pmax(diag(Jac %*% V_naive %*% t(Jac)), 0)), fit$K, fit$J)
  expect_true(mean(se_naive) > 3 * mean(m$se))         # sandwich << naive Hessian

  # reproduce the sandwich SE with an independent recomputation
  se_ref <- matrix(sqrt(pmax(diag(Jac %*% infometrics:::.mixed_gce_sandwich(fit) %*% t(Jac)), 0)),
                   fit$K, fit$J)
  expect_equal(unname(m$se), unname(se_ref), tolerance = 1e-6)
})

test_that("margins(se = FALSE) unchanged; se needs average = TRUE", {
  d   <- sim_mixed(N = 25L)
  fit <- mixed_gce(d$Y, d$X, nu = d$nu)
  expect_true(is.matrix(margins(fit)))
  expect_equal(dim(margins(fit, average = FALSE)), c(fit$N, fit$J, fit$K))
  expect_error(margins(fit, se = TRUE, average = FALSE), "average = TRUE")
})

test_that("margins bootstrap SEs run and roughly match the sandwich", {
  d   <- sim_mixed(N = 30L)
  fit <- mixed_gce(d$Y, d$X, nu = d$nu)
  set.seed(21)
  ms  <- margins(fit, se = TRUE, se_method = "sandwich")
  mb  <- margins(fit, se = TRUE, se_method = "bootstrap", B = 30L)
  expect_equal(mb$se_method, "bootstrap")
  expect_equal(dim(mb$se), c(fit$K, fit$J))
  expect_true(all(is.finite(mb$se) & mb$se >= 0))
  expect_true(max(ms$se / mb$se) < 3 && min(ms$se / mb$se) > 1 / 3)  # same order
})
