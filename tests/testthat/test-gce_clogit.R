# ============================================================
# test-gce_clogit.R
# Tests for gce_clogit() — conditional/mixed logit GCE estimator.
# ============================================================

# ---- Shared data simulator --------------------------------------------------
# Generates long-format Fishing-style data:
#   N individuals, J alternatives, price & catch (choice-specific),
#   income (individual-specific), choices via a simple conditional logit.

.sim_clogit <- function(N = 60L, seed = 1L) {
  set.seed(seed)
  alts <- c("beach", "boat", "charter", "pier")
  J    <- length(alts)
  df   <- data.frame(
    id     = rep(seq_len(N), each = J),
    alt    = rep(alts, times = N),
    price  = round(runif(N * J, 50, 300)),
    catch  = round(runif(N * J, 0.1, 1.5), 2),
    income = rep(round(runif(N, 20, 80)), each = J)
  )
  # True utility: beta_price = -0.01, beta_catch = 0.5,
  #              gamma_income_boat = 0.02 (example)
  eta   <- with(df, -0.01 * price + 0.5 * catch)
  probs <- tapply(exp(eta), df$id, function(x) x / sum(x))
  chosen <- unlist(lapply(seq_len(N), function(i)
    seq_len(J) == sample(J, 1L, prob = probs[[i]])))
  df$chosen <- chosen
  df
}


# ---- Basic convergence and structure ----------------------------------------

test_that("gce_clogit() converges on simulated data", {
  df  <- .sim_clogit(N = 60L, seed = 1L)
  fit <- gce_clogit(chosen ~ price + catch | income,
                    data = df, alt = "alt", id = "id")
  expect_true(fit$converged)
})

test_that("gce_clogit() produces K = (J-1) + K_cs + (J-1)*K_is parameters", {
  df  <- .sim_clogit(N = 60L, seed = 2L)
  fit <- gce_clogit(chosen ~ price + catch | income,
                    data = df, alt = "alt", id = "id")
  J    <- 4L   # beach, boat, charter, pier
  K_cs <- 2L   # price, catch
  K_is <- 1L   # income
  expect_equal(fit$K, (J - 1L) + K_cs + (J - 1L) * K_is)
  expect_equal(fit$J, J)
  expect_equal(length(fit$beta), fit$K)
})

test_that("gce_clogit() moment constraint: gradient ~ 0 at optimum", {
  df  <- .sim_clogit(N = 60L, seed = 3L)
  fit <- gce_clogit(chosen ~ price + catch | income,
                    data = df, alt = "alt", id = "id")
  # Gradient of the dual = sum_{ij}(Y-P-E)*Z, should be ~0
  grad_norm <- max(abs(.clogit_grad(
    fit$lambda, fit$y_mat, fit$Z, fit$v,
    log(fit$P0), log(fit$W0)
  )))
  expect_lt(grad_norm, 1e-3)
})


# ---- Formula variants -------------------------------------------------------

test_that("gce_clogit(): choice-specific only (no individual vars)", {
  df  <- .sim_clogit(N = 60L, seed = 4L)
  fit <- gce_clogit(chosen ~ price + catch | 0,
                    data = df, alt = "alt", id = "id")
  J    <- 4L; K_cs <- 2L
  expect_equal(fit$K, (J - 1L) + K_cs)
  expect_true(fit$converged)
})

test_that("gce_clogit(): individual-specific only (no choice vars)", {
  df  <- .sim_clogit(N = 60L, seed = 5L)
  fit <- gce_clogit(chosen ~ 0 | income,
                    data = df, alt = "alt", id = "id")
  J    <- 4L; K_is <- 1L
  expect_equal(fit$K, (J - 1L) + (J - 1L) * K_is)
  expect_true(fit$converged)
})

test_that("gce_clogit(): ASC-only model (intercepts only)", {
  df  <- .sim_clogit(N = 60L, seed = 6L)
  fit <- gce_clogit(chosen ~ 0 | 0,
                    data = df, alt = "alt", id = "id")
  expect_equal(fit$K, 3L)    # J - 1 = 3 ASCs
  expect_true(fit$converged)
})


# ---- Probability validity ---------------------------------------------------

test_that("gce_clogit() p_hat rows are valid probability vectors", {
  df  <- .sim_clogit(N = 60L, seed = 7L)
  fit <- gce_clogit(chosen ~ price + catch | income,
                    data = df, alt = "alt", id = "id")
  expect_lt(max(abs(rowSums(fit$p_hat) - 1)), 1e-8)
  expect_true(all(fit$p_hat > 0))
})

test_that("gce_clogit() w_hat (i,j) fibers are valid probability vectors", {
  df  <- .sim_clogit(N = 60L, seed = 8L)
  fit <- gce_clogit(chosen ~ price + catch | income,
                    data = df, alt = "alt", id = "id")
  w_sums <- apply(fit$w_hat, c(1, 2), sum)
  expect_lt(max(abs(w_sums - 1)), 1e-8)
  expect_true(all(fit$w_hat > 0))
})


# ---- Entropy bounds ---------------------------------------------------------

test_that("gce_clogit() S_p and S_w are in [0, 1]", {
  df  <- .sim_clogit(N = 60L, seed = 9L)
  fit <- gce_clogit(chosen ~ price + catch | income,
                    data = df, alt = "alt", id = "id")
  expect_gte(fit$S_p, 0); expect_lte(fit$S_p, 1)
  expect_gte(fit$S_w, 0); expect_lte(fit$S_w, 1)
})


# ---- GME == GCE with uniform priors -----------------------------------------

test_that("GME-clogit == GCE-clogit with explicit uniform priors", {
  df      <- .sim_clogit(N = 60L, seed = 10L)
  fit_gme <- gce_clogit(chosen ~ price + catch | income,
                        data = df, alt = "alt", id = "id")
  N <- fit_gme$N; J <- fit_gme$J; M <- fit_gme$M
  p0_unif <- matrix(1 / J, N, J)
  w0_unif <- array(1 / M, c(N, J, M))
  fit_gce <- gce_clogit(chosen ~ price + catch | income,
                        data = df, alt = "alt", id = "id",
                        p0 = p0_unif, w0 = w0_unif,
                        init = fit_gme$lambda)
  expect_lt(max(abs(coef(fit_gme) - coef(fit_gce))), 1e-6)
})


# ---- S3 methods -------------------------------------------------------------

test_that("gce_clogit() S3 methods return correct types", {
  df  <- .sim_clogit(N = 60L, seed = 11L)
  fit <- gce_clogit(chosen ~ price + catch | income,
                    data = df, alt = "alt", id = "id")
  expect_equal(length(coef(fit)),       fit$K)
  expect_equal(dim(fitted(fit)),        c(fit$N, fit$J))
  expect_equal(dim(residuals(fit)),     c(fit$N, fit$J))
  expect_output(print(fit))
  expect_output(summary(fit))
})

test_that("gce_clogit() coefficient names are correctly structured", {
  df  <- .sim_clogit(N = 60L, seed = 12L)
  fit <- gce_clogit(chosen ~ price + catch | income,
                    data = df, alt = "alt", id = "id")
  nms <- names(fit$beta)
  # ASCs: "(Intercept):boat", "(Intercept):charter", "(Intercept):pier"
  expect_true(any(grepl("^\\(Intercept\\):", nms)))
  # Generic coefs: "price", "catch"
  expect_true("price" %in% nms)
  expect_true("catch" %in% nms)
  # Alt-specific coefs: "income:boat", "income:charter", "income:pier"
  expect_true(any(grepl("^income:", nms)))
})

test_that("gce_clogit() beta = -lambda", {
  df  <- .sim_clogit(N = 60L, seed = 13L)
  fit <- gce_clogit(chosen ~ price + catch | income,
                    data = df, alt = "alt", id = "id")
  expect_equal(fit$beta, -fit$lambda)
})


# ---- Input validation -------------------------------------------------------

test_that("gce_clogit() error: p0 rows not summing to 1", {
  df     <- .sim_clogit(N = 60L, seed = 14L)
  p0_bad <- matrix(0.3, nrow = 60L, ncol = 4L)  # rows sum to 1.2
  expect_error(
    gce_clogit(chosen ~ price + catch | income,
               data = df, alt = "alt", id = "id", p0 = p0_bad),
    "Rows of p0 must sum to 1"
  )
})

test_that("gce_clogit() error: wrong number of rows in data", {
  df_bad <- .sim_clogit(N = 60L, seed = 15L)
  df_bad <- df_bad[-1L, ]   # drop one row -> unbalanced panel
  expect_error(
    gce_clogit(chosen ~ price + catch | income,
               data = df_bad, alt = "alt", id = "id"),
    "one row per"
  )
})

test_that("gce_clogit() error: wrong init length", {
  df  <- .sim_clogit(N = 60L, seed = 16L)
  expect_error(
    gce_clogit(chosen ~ price + catch | income,
               data = df, alt = "alt", id = "id",
               init = rep(0, 3)),
    "init must have length K"
  )
})


# ---- warm_start = FALSE gives same optimum ----------------------------------

test_that("warm_start = FALSE converges to same point as TRUE", {
  df     <- .sim_clogit(N = 60L, seed = 17L)
  fit_ws  <- gce_clogit(chosen ~ price + catch | income,
                        data = df, alt = "alt", id = "id",
                        warm_start = TRUE)
  fit_nws <- gce_clogit(chosen ~ price + catch | income,
                        data = df, alt = "alt", id = "id",
                        warm_start = FALSE)
  expect_true(fit_ws$converged)
  expect_true(fit_nws$converged)
  # Convex objective => same global minimum
  expect_lt(max(abs(coef(fit_ws) - coef(fit_nws))), 1e-4)
})
