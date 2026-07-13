# Tests for markov_ce() — Cross-Entropy Markov transition matrix from a panel
# (Golan 2008, Sec. 7.7.1).

# Vectorized one-hot panel simulator from a known transition matrix.
sim_markov_panel <- function(P, N, Tn, seed, init = NULL) {
  set.seed(seed)
  K <- nrow(P)
  s <- if (is.null(init)) sample(K, N, replace = TRUE) else init
  id <- rep(seq_len(N), each = Tn); tm <- rep(seq_len(Tn), N)
  state <- integer(N * Tn)
  for (t in seq_len(Tn)) {
    state[tm == t] <- s
    s <- vapply(s, function(k) sample(K, 1L, prob = P[k, ]), integer(1))
  }
  oh <- diag(K)[state, ]
  d <- data.frame(i = id, t = tm)
  for (k in seq_len(K)) d[[paste0("y", k)]] <- oh[, k]
  d
}

P_true <- matrix(c(.6, .3, .1, .2, .5, .3, .1, .3, .6), 3, byrow = TRUE)
st3 <- c("y1", "y2", "y3")

test_that("recovers a known transition matrix (indicator panel)", {
  panel <- sim_markov_panel(P_true, N = 250L, Tn = 6L, seed = 1)
  fit <- markov_ce(panel, id = "i", time = "t", states = st3)

  expect_s3_class(fit, "markov_ce")
  expect_s3_class(fit, "infometrics")
  expect_equal(fit$method, "dual")
  expect_equal(fit$state_type, "indicator")
  expect_lt(max(abs(coef(fit) - P_true)), 0.1)
  expect_lt(max(abs(rowSums(coef(fit)) - 1)), 1e-8)
  expect_lt(fit$moment_residual, 1e-4)
  expect_true(fit$converged)
})

test_that("CE with a uniform prior equals ME (p0 = NULL) on an identified panel", {
  panel <- sim_markov_panel(P_true, N = 200L, Tn = 5L, seed = 2)
  f_me <- markov_ce(panel, "i", "t", st3)
  f_ce <- markov_ce(panel, "i", "t", st3, p0 = matrix(1 / 3, 3, 3))
  expect_equal(f_me$p, f_ce$p, tolerance = 1e-8)
})

test_that("under-identified rows fall back to the prior", {
  # 2-state chain among {1,2}: state 3 is never an origin -> rank(A) < K.
  Q <- matrix(c(.6, .4, 0, .3, .7, 0, 0, 0, 0), 3, byrow = TRUE)
  panel <- sim_markov_panel(Q, N = 200L, Tn = 5L, seed = 5,
                            init = sample(1:2, 200L, replace = TRUE))

  expect_warning(g1 <- markov_ce(panel, "i", "t", st3), "rank")
  pr <- matrix(c(.5, .3, .2, .2, .5, .3, .15, .25, .60), 3, byrow = TRUE)
  g2 <- suppressWarnings(markov_ce(panel, "i", "t", st3, p0 = pr))

  expect_equal(g1$n_from[3], 0)                          # state 3 never an origin
  expect_equal(unname(g1$p[3, ]), rep(1 / 3, 3), tolerance = 1e-4)  # uniform fallback
  expect_equal(unname(g2$p[3, ]), pr[3, ], tolerance = 1e-4)        # prior fallback
  expect_lt(max(abs(g1$p[1:2, ] - g2$p[1:2, ])), 1e-4)  # identified rows unchanged
})

test_that("compositional shares are handled", {
  set.seed(3)
  N <- 150L; Tn <- 5L; K <- 3L
  sh <- matrix(runif(N * Tn * K), ncol = K); sh <- sh / rowSums(sh)
  panel <- data.frame(i = rep(seq_len(N), each = Tn), t = rep(seq_len(Tn), N),
                      y1 = sh[, 1], y2 = sh[, 2], y3 = sh[, 3])
  fit <- markov_ce(panel, "i", "t", st3)

  expect_equal(fit$state_type, "shares")
  expect_lt(max(abs(rowSums(coef(fit)) - 1)), 1e-8)
  expect_lt(fit$moment_residual, 1e-4)
})

test_that("covariates are accepted and the matrix stays row-stochastic", {
  panel <- sim_markov_panel(P_true, N = 200L, Tn = 5L, seed = 4)
  set.seed(9); panel$z <- rnorm(nrow(panel))
  fit <- suppressWarnings(markov_ce(panel, "i", "t", st3, covariates = "z"))
  expect_lt(max(abs(rowSums(coef(fit)) - 1)), 1e-6)
  expect_equal(dim(fit$lambda), c(1L, 3L))             # S x K
})

test_that("input validation throws informative errors", {
  panel <- sim_markov_panel(P_true, N = 80L, Tn = 4L, seed = 6)

  expect_error(markov_ce(panel[-1, ], "i", "t", st3), "balanced")
  bad <- panel; bad$y1[1] <- 0.5
  expect_error(markov_ce(bad, "i", "t", st3), "sum to 1")
  expect_error(markov_ce(panel, "i", "t", "y1"), "at least two state columns")
  expect_error(markov_ce(panel, "i", "t", st3, p0 = matrix(1 / 2, 2, 2)),
               "p0 must be K x K")
  expect_error(markov_ce(panel, "i", "t", c("y1", "y2", "nope")), "not found")

  one_period <- panel[panel$t == 1, ]
  expect_error(markov_ce(one_period, "i", "t", st3), "two periods")
})

test_that("S3 methods run and return correctly shaped objects", {
  panel <- sim_markov_panel(P_true, N = 150L, Tn = 5L, seed = 7)
  fit <- markov_ce(panel, "i", "t", st3)
  nt <- fit$n_transitions

  expect_equal(dim(coef(fit)), c(3L, 3L))               # K x K
  expect_equal(dim(fitted(fit)), c(nt, 3L))
  expect_equal(dim(residuals(fit)), c(nt, 3L))
  expect_output(print(fit), "Markov")
  s <- summary(fit)
  expect_s3_class(s, "summary.markov_ce")
  expect_output(print(s), "Transition matrix")
})

# ---- fano_bounds() and margins() (Golan 2008, sec 7.5 / 7.7) -----------------

Hbin <- function(pe) ifelse(pe > 0 & pe < 1, -pe * log(pe) - (1 - pe) * log(1 - pe), 0)

test_that("fano_bounds returns per-origin bounds that hold", {
  panel <- sim_markov_panel(P_true, N = 250L, Tn = 6L, seed = 1)
  fit <- markov_ce(panel, "i", "t", st3)
  fb <- fano_bounds(fit)
  expect_s3_class(fb, "data.frame")
  expect_equal(nrow(fb), length(st3))                 # one row per origin state
  expect_named(fb, c("p_max", "pe", "H", "S", "pe_lower"))
  expect_equal(fb$pe, unname(1 - apply(fit$p_hat, 1, max)), tolerance = 1e-12)
  expect_true(all(fb$S >= 0 & fb$S <= 1 + 1e-12))
  expect_true(all(fb$pe >= fb$pe_lower - 1e-9))                            # weak
  expect_true(all(Hbin(fb$pe) + fb$pe * log(length(st3) - 1) >= fb$H - 1e-9))  # strong
  expect_named(attr(fb, "overall"), c("mean_pe", "mean_pe_lower", "S_system"))
})

test_that("margins requires covariates; no-covariate fit errors", {
  panel <- sim_markov_panel(P_true, N = 200L, Tn = 5L, seed = 2)
  fit <- markov_ce(panel, "i", "t", st3)
  expect_null(fit$se_data)
  expect_error(margins(fit), "require covariates")
})

test_that("margins matches the softmax-Jacobian identity (covariate fit)", {
  panel <- sim_markov_panel(P_true, N = 200L, Tn = 6L, seed = 3)
  set.seed(42)
  for (s in 1:3) panel[[paste0("z", s)]] <- rnorm(nrow(panel)) + panel[[paste0("y", s)]]
  fit <- suppressWarnings(markov_ce(panel, "i", "t", st3, covariates = c("z1", "z2", "z3")))
  ME <- margins(fit, average = FALSE)                 # K x K x S
  expect_equal(dim(ME), c(3L, 3L, 3L))
  P <- fit$p_hat; lam <- fit$lambda
  for (k in 1:3) { Jk <- diag(P[k, ]) - outer(P[k, ], P[k, ])
    for (s in 1:3) expect_equal(unname(ME[k, , s]), unname(as.vector(Jk %*% lam[s, ])),
                                tolerance = 1e-10) }
  expect_equal(dim(margins(fit)), c(3L, 3L))          # AME S x K
})

# ---- analytic-Hessian standard errors of lambda (Golan 2008 eq. 4.7) ---------

test_that("analytic Hessian equals the finite-difference Hessian", {
  panel <- sim_markov_panel(P_true, N = 200L, Tn = 6L, seed = 5)
  set.seed(9)
  for (s in 1:3) panel[[paste0("z", s)]] <- rnorm(nrow(panel)) + panel[[paste0("y", s)]]
  fit <- suppressWarnings(markov_ce(panel, "i", "t", st3, covariates = c("z1", "z2", "z3")))
  sd <- fit$se_data
  Han <- infometrics:::.markov_hessian(sd$A, fit$p_hat)
  grad <- function(lvec) as.vector(sd$M - crossprod(sd$A, infometrics:::.markov_P_eps(lvec, sd)$P))
  Hfd <- -infometrics:::.fd_jac(grad, as.vector(fit$lambda * sd$cov_scale))
  expect_lt(max(abs(Han - Hfd)), 1e-6)
})

test_that("se_lambda: shape, reference-column NA, and finite (indicator fit)", {
  panel <- sim_markov_panel(P_true, N = 250L, Tn = 6L, seed = 1)
  fit <- markov_ce(panel, "i", "t", st3)
  expect_equal(dim(fit$se_lambda), c(3L, 3L))
  expect_equal(fit$lambda_ref, 1L)
  expect_true(all(is.na(fit$se_lambda[, fit$lambda_ref])))        # reference column NA
  expect_true(all(is.finite(fit$se_lambda[, -fit$lambda_ref])))   # others finite
  expect_true(all(fit$se_lambda[, -fit$lambda_ref] >= 0))
  expect_output(print(summary(fit)), "analytic-Hessian sandwich")
})

test_that("analytic se_lambda matches a unit bootstrap (sandwich, not naive H^-1)", {
  panel <- sim_markov_panel(P_true, N = 250L, Tn = 6L, seed = 1)
  fit <- markov_ce(panel, "i", "t", st3); ref <- fit$lambda_ref
  units <- unique(panel$i); B <- 120L; LB <- array(NA_real_, c(3L, 3L, B))
  set.seed(3)
  for (b in seq_len(B)) {
    bu <- sample(units, length(units), replace = TRUE)
    dl <- do.call(rbind, lapply(seq_along(bu), function(m) { d <- panel[panel$i == bu[m], ]; d$i <- m; d }))
    fb <- tryCatch(markov_ce(dl, "i", "t", st3), error = function(e) NULL)
    if (!is.null(fb) && fb$converged) LB[, , b] <- fb$lambda - fb$lambda[, ref]
  }
  se_boot <- apply(LB, c(1, 2), function(z) sd(z, na.rm = TRUE))
  r <- fit$se_lambda[, -ref] / se_boot[, -ref]
  expect_true(median(r) > 0.6 && median(r) < 1.6)                 # sandwich ~ bootstrap

  # the naive information-matrix SE (H^-1) grossly overstates (~sqrt(n_from))
  sd <- infometrics:::.markov_hessian(diag(fit$n_from), fit$p_hat)
  free <- which(rep(1:3, each = 3) != ref)
  se_naive <- matrix(NA_real_, 3, 3)
  se_naive[, -ref] <- matrix(sqrt(diag(solve(sd[free, free]))), 3, 2)
  expect_gt(median(se_naive[, -ref] / fit$se_lambda[, -ref]), 5)  # H^-1 >> sandwich
})

test_that("a (near-)boundary / non-converged covariate fit yields NA se_lambda without error", {
  panel <- sim_markov_panel(P_true, N = 400L, Tn = 5L, seed = 3)
  set.seed(1)
  for (s in 1:3) panel[[paste0("z", s)]] <- rnorm(nrow(panel)) + 3 * panel[[paste0("y", s)]]
  fit <- suppressWarnings(markov_ce(panel, "i", "t", st3, covariates = c("z1", "z2", "z3")))
  expect_true(is.matrix(fit$se_lambda))                          # no error
  expect_equal(dim(fit$se_lambda), c(3L, 3L))
})
