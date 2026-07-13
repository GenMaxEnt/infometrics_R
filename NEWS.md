# infometrics 0.3.0

## Package scope reduced for CRAN submission

The package now focuses on the Generalized Maximum Entropy (GME) and
Generalized Cross-Entropy (GCE) estimator family. The following functions were
**removed**: `me()`, `gme()`, `me_mnl()`, `gme_mnl()`, `gce_clogit()`,
`gce_clogitWrap()`, and `gce_table()`. The retained estimators are
`inverse_ce()`, `inverse_noise()`, `linreg()`, `linreg_iv()`, `panel_gce()`,
`matrix_ce()`, `matrix_gce()`, `markov_ce()`, `markov_gce()`,
`multinomial_gce()`, and `mixed_gce()`.

The `entropy.R` module was trimmed to `shannon_entropy()`; the divergence
measures `kl_divergence()`, `renyi_entropy()`, `renyi_divergence()`,
`tsallis_divergence()`, `cressie_read()`, and `normalized_entropy()` were
removed. The deprecated `inverse_pure()` alias was also dropped &mdash; use
`inverse_ce()` directly.

### `inverse_ce()` (was `inverse_pure()`)

* `inverse_pure()` has been **renamed to `inverse_ce()`** (the method is
  cross-entropy: a uniform prior gives ME, a non-uniform prior gives CE).
* `inverse_ce()` now reports **information-matrix standard errors**
  `Var(lambda) = I^-1` (`se_lambda`, `vcov()`) and delta-method `se_p`, plus
  **`fano_bounds()`** for the recovered `p`. These are curvature/identification
  quantities (Golan 2008, sec 4.2), not sampling SEs; `NA`/`NULL` when the
  information matrix is singular.
* The prior argument was **renamed from `q` to `p0`** (the GME/GCE signal-prior
  name; for `inverse_ce` it is a K-vector prior over the states).
* Output cleanup: removed duplicate fields (kept `p_hat`/`lambda_hat`/`objective`/
  `p0`); dropped `converged` (the raw `convergence` code is kept).

### `inverse_noise()`

* **Standard errors** for `lambda_hat` and `p_hat` via `se_method` /
  `vcov(type=)` / `summary(se_method=)`: `"sandwich"` (default, robust HC0 —
  accurate), `"delta"` (naive `H^-1`; overstates ~10x, for comparison only), and
  `"bootstrap"` (residual bootstrap; stable aligned SEs for both parameters).
  **`vcov()`'s default changed** from the raw `H^-1` to the sandwich, since
  `H^-1` is not a valid sampling covariance here (`type = "delta"` recovers it).
* `summary()` now prints a `p` coefficient table with **standard errors and
  t-stats**, reports **normalized entropy for both signal `S(p)` and noise
  `S(w)`**, and adds a **Fano** line; **`fano_bounds()`** is available for
  `p_hat`.
* Output cleanup: removed duplicate fields (kept `p_hat`/`lambda_hat`/`w_hat`/
  `residuals`/`objective`); dropped `converged` (the raw `convergence` code is
  kept).

# infometrics 0.2.0

## Phase 2 — GME/GCE regression estimator

### New functions

**GME / GCE estimator** (`R/gme.R`):

* `gme(formula, data, Z, p0, V, w0, nu, method, control)` — Generalized
  Maximum Entropy (GME) and Generalized Cross-Entropy (GCE) estimation for
  linear regression, following Golan, Judge and Miller (1996) and Golan
  (2008, Chapter 6). Uniform priors (`p0 = NULL`, `w0 = NULL`) give GME;
  user-supplied priors give GCE. One function, not two — the same relationship
  as `me()` to ME/CE.
* The dual concentrated problem (Golan 2008, below Eq. 6.9) is minimised via
  BFGS with analytic gradient. Both signal (`Omega_k`) and error (`Psi_t`)
  partition functions use the log-sum-exp trick for numerical stability.
* S3 methods: `print`, `summary`, `coef`, `fitted`, `residuals`, `vcov`.
* Returns per-coefficient normalised entropy `S_p`, overall `S_P`, and
  `pseudo_R2` (Golan 2008, Section 7.5).

### References

Golan, A., Judge, G. and Miller, D. (1996). *Maximum Entropy Econometrics*.
Wiley.

---

# infometrics 0.1.0

## Phase 1 — Initial release

### New functions

**Entropy measures** (`R/entropy.R`):

* `shannon_entropy(p, base)` — Shannon entropy H(p) = -Σ p log p, with
  configurable base (nats, bits, hartleys). Convention 0·log(0) = 0 applied.
* `kl_divergence(p, q, base)` — Kullback-Leibler divergence D(p‖q).
* `renyi_entropy(p, alpha, base)` — Rényi entropy of order α; reduces to
  Shannon as α → 1.
* `renyi_divergence(p, q, alpha, base)` — Rényi cross-entropy of order α
  (Golan 2008, Eq. 3.6).
* `tsallis_divergence(p, q, alpha)` — Tsallis cross-entropy of order α
  (Golan 2008, Eq. 3.7).
* `cressie_read(p, q, alpha)` — Cressie-Read power divergence, the unifying
  criterion for the IT estimator class (Golan 2008, Eq. 3.8). Special cases
  include KL (α → 0), Pearson χ² (α = 1), EL (α → -1).
* `normalized_entropy(p, q)` — S(p̃) = H(p) / H(q), the primary
  goodness-of-fit statistic for IT models (Golan 2008, Section 6.4).

**ME / CE estimator** (`R/me.R`):

* `me(y, X, q, method, control)` — Maximum Entropy (ME) and Cross-Entropy (CE)
  estimation for under-determined systems. Both primal (constrained, via
  `nloptr`) and dual (concentrated, unconstrained, default) formulations
  are implemented.
* S3 methods: `print`, `summary`, `coef`, `fitted`, `residuals`.
* The dual (concentrated) form solves the T-dimensional unconstrained problem
  over Lagrange multipliers λ via BFGS, reducing dimensionality from K >> 1
  to T. No extra package dependency needed.

**Utilities** (`R/utils.R`):

* `normalize_data(x, by)` — rescales data to [0, 1] to prevent overflow in
  partition function computations.
* `make_support(half_range, M, center)` — constructs symmetric M-point support
  spaces for GME/GCE estimation.
* `default_supports(y, X, ...)` — constructs data-driven default signal (Z)
  and error (V) support spaces from OLS estimates, using the three-sigma rule
  for error bounds (Pukelsheim 1994).

### References

Golan, A. (2008). Information and Entropy Econometrics — A Review and
Synthesis. *Foundations and Trends in Econometrics*, **2**(1-2), 1-145.
