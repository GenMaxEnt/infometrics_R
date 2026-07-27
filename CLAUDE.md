# CLAUDE.md — infometrics

This file instructs Claude on how to work within the `infometrics` R package
project. Read it fully before making any changes.

---

## Project overview

`infometrics` is an R package implementing the class of **Information-Theoretic
(IT) estimators** for econometric models, following the unified framework of
Golan (2008) *Information and Entropy Econometrics — A Review and Synthesis*
(Foundations and Trends in Econometrics, 2(1-2), 1–145). The package is
being built for submission to CRAN.

The central unifying concept is the **Cressie-Read power divergence criterion**
(Golan 2008, Eq. 3.8–3.9), a generalized entropy of order α that nests
Shannon entropy (α→0), Empirical Likelihood (α→−1), and the Pearson
chi-squared statistic (α=1) as special cases.

---

## Current state

### Phase 1 — COMPLETE (`v0.1.0`)
All Phase 1 code lives in `R/`, tests in `tests/testthat/`.

| File | Contents | Status |
|---|---|---|
| `R/entropy.R` | `shannon_entropy` (+ internal `.check_prob`); trimmed to the one measure the retained estimators depend on | ✅ complete |
| `R/inverse_ce.R` | `inverse_ce()` — ME/CE pure inverse problem via the dual concentrated model, **formula interface** (lm-style); now reports **information-matrix SEs** `Var(λ)=I⁻¹` (`se_lambda`, `vcov`) + delta-method `se_p` and **`fano_bounds`** for the recovered `p`; S3 class `c("inverse_ce", "infometrics")` with `print`, `summary`, `coef`, `fitted`, `residuals`, `vcov`, `fano_bounds` (+ `print.summary.inverse_ce`). De-duplicated output (keeps `p_hat`/`lambda_hat`/`objective`/`p0`; dropped `converged`); prior arg is **`p0`** (renamed from `q`). `inverse_pure()` is a **deprecated alias** (still accepts legacy `q`). See "inverse_ce" below | ✅ complete, 68/68 tests passing |
| `R/matrix_ce.R` | `matrix_ce()` — CE/ME **matrix balancing** (Golan §7.2): recovers an n×m column-stochastic matrix from column weights `x` and row aggregates `y` (`y = Px`, `colSums = 1`) via the dual; S3 class `c("matrix_ce", "infometrics")` with `print`, `summary` (Golan §7.5 information measures + ER test), `coef`, `fitted`, `residuals` (+ `print.summary.matrix_ce`). **Genuinely new estimator** (not a duplicate) — see "matrix_ce" below | ✅ complete, 42/42 tests passing |
| `R/markov_ce.R` | `markov_ce()` — CE estimation of a first-order **Markov transition matrix** from a balanced long panel (Golan §7.7.1): K×K row-stochastic `P` from one-hot or compositional-share states, optional lead/lagged covariates; internal `.markov_arrays()` panel reshaper; S3 class `c("markov_ce", "infometrics")` with `print`, `summary`, `coef`, `fitted`, `residuals` (+ `print.summary.markov_ce`). **Genuinely new estimator** — see "markov_ce" below | ✅ complete, 31/31 tests passing |
| `R/utils.R` | `normalize_data`, `make_support`, `default_supports` | ✅ complete |
| `DESCRIPTION` | Package metadata, v0.3.0 | ✅ complete |
| `NAMESPACE` | Exports and S3 registrations | ✅ complete |
| `NEWS.md` | Changelog | ✅ complete |

### Phase 2 — COMPLETE (`v0.2.0`)

| File | Contents | Status |
|---|---|---|
| `R/inverse_noise.R` | `inverse_noise()` — noisy (stochastic-moment) GME/GCE inverse problem via the dual concentrated model, **formula interface** (lm-style); now reports **sampling SEs** for `lambda_hat`/`p_hat` by `se_method`/`vcov(type=)`/`summary(se_method=)` ∈ {`sandwich` (default), `delta`, `bootstrap`}, **signal + noise normalized entropy** (`S`/`S_w`), and **`fano_bounds`** for `p_hat`; de-duplicated output (keeps `p_hat`/`lambda_hat`/`w_hat`/`residuals`/`objective`; dropped `converged`); S3 class `c("inverse_noise", "infometrics")` with `print`, `summary` (p-table w/ SE+t), `coef`, `fitted`, `residuals`, `vcov`, `fano_bounds` (+ `print.summary.inverse_noise`). See "inverse_noise" below | ✅ complete, 80/80 tests passing |
| `tests/testthat/test-inverse_noise.R` | 80 tests covering convergence, noisy moment condition (residuals = recovered noise), shrink-to-pure limit, GME vs GCE, support-width effect, PD Hessian, de-dup (dropped aliases absent), signal+noise normalized entropy, SEs (sandwich stored/finite, vcov(type=) shapes, delta≫sandwich, bootstrap≈sandwich band, se_method=none), Fano bounds (weak+strong), scalar-v support, input validation, S3 methods (p-table SE+t, entropy+Fano lines, alt se_method) | ✅ complete |
| `R/linreg.R` | `linreg()` — GCE **linear regression** (`y = Xβ + e`), **formula interface** (lm-style) with `predict`, `r.squared`, and an **entropy-ratio (ER) `summary`** (per-coefficient χ²(1) ER test + overall slopes test, via restricted refits); S3 class `c("linreg", "infometrics")` with `coef`, `fitted`, `residuals`, `vcov`, `predict`, `print`, `summary` (+ `print.summary.linreg`). Internal solver factored into `.linreg_engine()` (reused by the ER refits) with `.linreg_Hstar()`. **Same estimator as `gme()`** — see "linreg vs gme" below | ✅ complete, 51/51 tests passing |
| `tests/testthat/test-linreg.R` | 51 tests covering convergence, OLS limit, tight-support shrinkage, moment condition, GME vs GCE, shared-vs-matrix support, predict, PD Hessian/vcov, entropy/R², ER-test validity, Golan §6.6 replication (collinearity), overall ER test, input validation, S3 methods | ✅ complete |
| `R/matrix_gce.R` | `matrix_gce()` — GCE **matrix balancing with stochastic moments** (Golan §7.4): noisy sibling of `matrix_ce()`. Recovers a column-stochastic signal `p` (n×m) **and** a row-stochastic noise `w` (n×H) from `y = Px + e`, `e_i = Σ_h v_h w_ih`, with `nu` weight; S3 class `c("matrix_gce", "infometrics")` with `coef`, `fitted`, `residuals`, `print`, `summary`, `print.summary`. **Genuinely new estimator** — see "matrix_gce" below | ✅ complete, 33/33 tests passing |
| `tests/testthat/test-matrix_gce.R` | 33 tests covering convergence, stochastic moment condition, matrix_ce shrink-limit, GME vs GCE, entropy bounds, §7.5 summary (W/df/p-value, per-column S), input validation, normalization warning, S3 methods | ✅ complete |
| `R/markov_gce.R` | `markov_gce()` — GCE (noisy-moment) **Markov transition matrix** from a panel with covariates (Golan §7.7.1, eq. 7.28): noisy sibling of `markov_ce()`. Per-moment error on a symmetric `[−1,1]` support `v`, weight `nu`; reuses `.markov_arrays()`; S3 class `c("markov_gce", "infometrics")` with `print`, `summary`, `coef`, `fitted`, `residuals` (+ `print.summary.markov_gce`). **Genuinely new estimator** — see "markov_gce" below | ✅ complete, 28/28 tests passing |
| `tests/testthat/test-markov_gce.R` | 28 tests covering noisy-FOC convergence, nu/v effects (toward prior), canonical aliases, input validation, S3 methods | ✅ complete |
| `R/multinomial_gce.R` | `multinomial_gce()` — **nu-weighted GCE multinomial** (Golan-Judge-Perloff 1996), **matrix `(y, X)` interface**; recovers N×J signal `p` + N×J×M noise `w`; `which_alternative` reference; new exported **`margins`** generic (∂p_ij/∂x_ik) **now with sandwich/delta/bootstrap SEs** (`margins(fit, se=TRUE)`, default robust sandwich) and a new exported **`fano_bounds`** generic — **Fano-inequality error bounds** for `p_hat` (Golan §7.5). S3 class `c("multinomial_gce", "infometrics")` with `coef` (→ `lambda`), `fitted`, `residuals`, `print`, `summary` (now w/ a Fano line), `margins`, `fano_bounds`. **Overlaps `gme_mnl`** (≡ it at nu=0.5/matched support/ref 1) — see "multinomial_gce vs gme_mnl" below | ✅ complete, 62/62 tests passing |
| `tests/testthat/test-multinomial_gce.R` | 62 tests covering convergence, exact `gme_mnl` reproduction, recovery, nu/v effects, prior shift, `which_alternative` near-invariance, input validation, S3 methods, canonical aliases, `margins` vs a numerical derivative, **Fano bounds (weak + strong inequalities hold), and margins SEs (sandwich default ≈ bootstrap, delta identity, sandwich < delta ordering)** | ✅ complete |
| `R/linreg_iv.R` | `linreg_iv()` — **stochastic-moments GME-IV** regression (Golan §pp.89-91): IV sibling of `linreg`/`gme`, moment `IV'(y−Xβ−e)=0`, `β = Zp`; supports over-identification; S3 class `c("linreg_iv", "infometrics")` with `coef` (→ β), `fitted`, `residuals`, `print`, `summary`. **Genuinely new estimator** — see "linreg_iv" below | ✅ complete, 29/29 tests passing |
| `tests/testthat/test-linreg_iv.R` | 29 tests covering convergence/moment satisfaction, OLS recovery (IV=X), 2SLS recovery (endogeneity correction), over-identification, support shrinkage, input validation, canonical aliases, S3 methods | ✅ complete |
| `R/panel_gce.R` | `panel_gce()` — dual **GME/GCE for the one-way error-components panel model** `y_nt = x_nt'β + μ_n + ε_nt` (Lee & Cheon 2014, eq. 3.12): the panel/individual-effects sibling of `linreg`. Three reparameterized blocks — coefficients `β = Zp` (support `Z`, prior `p0`), **individual effects** `μ_n = Σ_r f_nr g_nr` (support `FF`, prior `g0`), errors `e_nt = Σ_j v_j w_ntj` (support `v`, prior `w0`) — solved via the concentrated dual over `λ` (length N·T), weight `nu`; `period`/`unit` layouts. **Standard errors** (Golan §3.3) for β/p/μ via `vcov_type` + delta method; S3 class `c("panel_gce", "infometrics")` with `coef` (→ β), `vcov`, `fitted`, `residuals`, `print`, `summary` (ER-style coef table w/ z & p-values). **Genuinely new estimator** — see "panel_gce" below | ✅ complete, 44/44 tests passing |
| `tests/testthat/test-panel_gce.R` | 44 tests covering FOC/convergence, within (fixed-effects) recovery (both layouts), individual-effects recovery (`cor(μ̂,μ)>0.9`), directional prior shift + narrow-support shrinkage, **SE surface (vcov PD/symmetric, se_beta/se_p/se_mu finite, μ finite-T floor, vcov_type ordering, se_p delta identity, p-value column)**, input validation, canonical aliases + S3 methods | ✅ complete |
| `R/mixed_gce.R` | `mixed_gce()` — **doubly-reparameterized GCE "mixed" model**: response `Y` (N×J shares), design `X` (N×J×K); the signal probability is *itself* reparameterized on a support `s∈[0,1]` (`θ_ijm ∝ θ0 exp(s_m(V_ij+ρ_i)/ν)`, `V_ij=Σ_k X_ijk λ_kj`, `p_ij=Σ_m s_m θ_ijm`), adding-up `Σ_j p_ij=1` via a per-obs multiplier `ρ_i`; noise `w_ijh ∝ w0 exp(u_h V_ij/(1−ν))`. Concentrated dual maximized over `(ρ,λ)`; flexible `s`/`u` (NULL/count/vector); dual-Hessian SEs. S3 class `c("mixed_gce", "infometrics")` with `coef` (→ λ), `fitted`, `residuals`, `vcov`, `print`, `summary`, `margins` (**now with robust-sandwich SEs**), `fano_bounds` (**Fano error bounds for p_hat**, Golan §7.5). **Genuinely new estimator** — see "mixed_gce" below | ✅ complete, 72/72 tests passing |
| `tests/testthat/test-mixed_gce.R` | 72 tests covering FOC/convergence, normalization + probability/fiber validity, moment satisfaction, shrink-to-pure recovery of λ*/Pstar, GME==GCE-uniform + prior shift, flexible s/u forms, dual-Hessian SEs (PD/finite), margins vs central-difference derivative, **Fano bounds (weak+strong), margins sandwich SEs (≪ naive Hessian, ≈ bootstrap)**, input validation, canonical aliases + S3 methods | ✅ complete |

Uniform priors (`p0 = NULL`, `w0 = NULL`) → GME (max entropy).
User-supplied priors → GCE (min cross-entropy). One function, not two.

**Scope note (v0.3.0).** The standalone `me()`, `gme()`, `me_mnl()`,
`gme_mnl()`, `gce_clogit()`, `gce_clogitWrap()`, and `gce_table()` functions were
**removed** for the CRAN submission, along with the divergence measures in
`entropy.R` (only `shannon_entropy` is retained) and the deprecated
`inverse_pure()` alias. The retained estimators are `inverse_ce`,
`inverse_noise`, `linreg`, `linreg_iv`, `panel_gce`, `matrix_ce`, `matrix_gce`,
`markov_ce`, `markov_gce`, `multinomial_gce`, and `mixed_gce`. Multinomial
GME/GCE is now provided by `multinomial_gce()` (matrix interface) and
`mixed_gce()`.

### Phase 3 — NOT STARTED
Target: `R/inference.R`
- Entropy-ratio test, Wald test, large deviation bounds
- `predict()` / `confint()` methods for the retained estimator classes
  (`linreg` already provides `predict`).

### Phase 4 — NOT STARTED
CRAN submission preparation: extended tests, `R CMD check` clean pass,
full vignette with GME replication examples.

---

## Architecture decisions (do not change without discussion)

### Dual (concentrated) model as default
All estimators solve the **dual (unconstrained) concentrated model** by
default, working in Lagrange multiplier space (dimension T) rather than
probability space (dimension K >> T). This follows Golan (2008, Section 4.2,
Eq. 4.4–4.5) and is computationally superior. All retained estimators use the
dual form exclusively.

The dual objective for ME/CE is:
```
min_λ  { -λ'y + log(Ω(λ)) }
where  Ω(λ) = Σ_k q_k * exp(Σ_t λ_t * x_{tk})
```

For GME/GCE (Phase 2), the dual is:
```
min_λ  { Σ_t y_t λ_t + Σ_k log Ω_k(λ) + Σ_t log Ψ_t(λ) }
where  Ω_k(λ) = Σ_m p0_km * exp(−z_km * Σ_t λ_t x_tk)   [signal partition fn]
       Ψ_t(λ) = Σ_j w0_tj * exp(−λ_t * v_j)              [error  partition fn]
```
When `p0_km = 1/M` and `w0_tj = 1/J` (uniform), this reduces to GME.
With user-supplied priors it is GCE. The dual form is identical in both
cases — only the partition functions change.

### Numerical stability
Always use the **log-sum-exp trick** when computing log(Ω(λ)):
```r
max_e <- max(exponents)
log_omega <- log(sum(q * exp(exponents - max_e))) + max_e
```
Never compute `exp(exponents)` directly without this stabilization.

### Optimizer choice
- **Dual model**: `stats::optim(method = "BFGS")` with analytic gradient.
  No extra package dependency. This is the solver for every estimator.

### S3 class hierarchy
All estimator output objects inherit from `"infometrics"`:
- `inverse_ce()` returns class `c("inverse_ce", "infometrics")`
- `inverse_noise()` returns class `c("inverse_noise", "infometrics")`
- `linreg()` returns class `c("linreg", "infometrics")`
- `linreg_iv()` returns class `c("linreg_iv", "infometrics")`
- `panel_gce()` returns class `c("panel_gce", "infometrics")`
- `matrix_ce()` returns class `c("matrix_ce", "infometrics")`
- `matrix_gce()` returns class `c("matrix_gce", "infometrics")`
- `markov_ce()` returns class `c("markov_ce", "infometrics")`
- `markov_gce()` returns class `c("markov_gce", "infometrics")`
- `multinomial_gce()` returns class `c("multinomial_gce", "infometrics")`
- `mixed_gce()` returns class `c("mixed_gce", "infometrics")`

Each estimator subsumes GCE via its signal/noise prior arguments (`p0`, `w0`):
uniform priors → GME (max entropy), user priors → GCE (min cross-entropy). One
function per estimator family — there is no separate `gce_*` class.

**`inverse_ce()`** — solves the pure-moment ME/CE dual (Golan §4.2) through an
`lm`-style formula interface: the response is the moment vector and each RHS
term is a state (use `- 1` to drop the intercept). *Warns* (does not error) when
`K ≤ T + 1`. Adds the analytic T×T dual Hessian `Cov_p(moments)` (Eq. 4.7). Its
prior arg is **`p0`** (a K-vector prior over the states — the same *role* as the
GME/GCE signal prior, though a vector rather than a K×M matrix). `inverse_ce`
minimises `-λ'y + log Ω(λ)` without `fnscale`. Uniform prior ⇒ ME, non-uniform
⇒ CE (hence the name).

**Information-matrix SEs + Fano bounds for `inverse_ce` (Golan §4.2/§7.5, from
`pure_inverse.Rmd`)** — two honest add-ons that use quantities the estimator
already has, both explicitly framed as **NOT sampling SEs** (a pure inverse
problem is deterministic). (1) **`Var(λ) = I⁻¹(λ)`** — the stored dual Hessian
`hessian` *is* the Fisher information `I(λ) = Cov_p(moments)` (Eq. 4.7; verified
== finite-diff), so the object now carries `vcov_lambda` (via `vcov()`),
`se_lambda = sqrt(diag)`, and delta-method `se_p` (`Var(p)=J I⁻¹ Jᵀ`,
`J_{k,t}=p_k(x_tk−E_p[x_t]) = ∂p_k/∂λ_t`). These are **information-matrix /
curvature** quantities measuring how well the moments identify λ (the Rmd's "λ_t
= marginal information content of moment t"; a moment with larger `Var_p(x_t)`
has a **smaller** `se_lambda_t` — verified). The inverse is **rank-checked via
`eigen`** (`.rank_checked_inverse`, base R, no MASS): when `I` is singular — a
**constant/collinear/affinely-dependent moment row** or `T ≥ K` (rank ≤ K−1) —
`vcov_lambda` is `NULL` and `se_lambda`/`se_p` are `NA` (honest, vs a
pseudo-inverse's misleadingly finite value). `summary()` prints a p-table
(`Estimate | Std. Error`) and a λ-table (`Estimate | Std. Error | z value`, no
p-value) with the "curvature not sampling" note; for sampling inference use
`inverse_noise()`. (2) **`fano_bounds.inverse_ce`** — the recovered `p`
is one distribution over K states, so Golan's §7.5 Fano bound applies directly
(the companion to `S`): modal error `pe = 1−max_k p_k ≥ S(p)−ln2/ln(K)` with
`S(p)=H(p)/ln(K)` (uniform-reference, **distinct** from `inverse_ce`'s
prior-relative reported `S=H(p)/H(p0)`). Reuses the shared `.fano_row_bounds` on a
1-row matrix; returns a one-row df; `summary()` prints a Fano line. See
[[gme-se-validation]].

**`matrix_ce()`** — cross-entropy **matrix balancing** (Golan §7.2), a
**genuinely new estimator** (not a duplicate of any other). Recovers an n×m
**column-stochastic** matrix `p` (each column a distribution over the n rows)
from column weights `x` (length m) and row-weighted aggregates `y` (length n):
`y_i = Σ_j p_ij x_j`, `Σ_i p_ij = 1`, solution `p̂_ij = p0_ij exp(λ_i x_j)/Ω_j`,
the dual minimised over the n row multipliers. Typically under-determined (n·m
unknowns, n + m constraints) → returns the min-cross-entropy (max-entropy if
`p0` uniform) matrix consistent with the data. **By user instruction it keeps
its own argument names** (`y`, `x`, `p0`) rather than `me()`'s `q`. Uniform `p0`
→ ME; user `p0` → CE. Minimises the convex dual `-Σ λ_i y_i + Σ_j log Ω_j` (no
`fnscale`), matching `me()`/`inverse_ce`. Carries `p_hat` (no `w_hat` — pure
CE, no noise term); `H_signal` is the scalar `H(p̂)`; normalized entropy is
`S = H(p̂)/H(p0)` (prior-relative, Golan §7.5; for uniform `p0` this is
`H/(m·log n) ∈ [0,1]`, `=1` when every column is uniform; for a non-uniform CE
prior it is relative to that prior). `summary()` reports the §7.5 information
measures — information index `I = 1−S`, pseudo-R² `= 1−S`, per-column
`S(p_j) = H(p_j)/H(p0_j)`, and the entropy-ratio test
`W = 2[H(p0)−H(p̂)] ~ χ²(n−1)` for H0: P = P0 — matching `matrix_gce` (object
also carries `entropy_p0 = H(p0)`). Note: when the data is exactly consistent
with uniform `p` (optimum `λ*=0`) `optim` can report `convergence=1` though `p̂`
is recovered exactly.

**`markov_ce()`** — CE estimation of a first-order, stationary **Markov
transition matrix** `P` (K×K, row-stochastic, `P[k,j] = k→j`) from a **balanced
long panel** (Golan §7.7.1). States are simplex-valued — one-hot membership
indicators or compositional shares; optional covariates enter lead/lagged. Long
data (one row per unit-period) is reshaped internally by `.markov_arrays()` to
N×Tn×K arrays; transitions `t−1→t` form cross moments `A = ΣX'Z2`, `M = ΣZ1'Y`;
the dual over the S×K multipliers gives `P̂_kj ∝ p0_kj exp(Σ_s A_ks λ_sj)`. **By
user instruction it keeps its own argument names** (`data`, `id`, `time`,
`states`, `covariates`, `p0`). **By explicit user choice the dual is kept in the
maximise form (`fnscale=-1`)** — the *one* estimator in the package that
maximises rather than minimises (validated identical to the minimise flip). The
loop variable was renamed `T → Tn` (avoids shadowing `TRUE`). Carries `p_hat`
(no `w_hat` — pure CE); `H_signal` is the scalar mean row entropy; `S =
H̄(p̂)/H̄(p0)` is prior-relative. **Key identification property:** when the
conditioning is full rank (`A = X'X`, the common indicator case), `A'P = M`
*exactly* determines `P` as the empirical/regression transition matrix, so the
prior is irrelevant for identified origin states; only **under-identified** rows
(`rank(A) < K`, e.g. a state never observed as an origin) fall back to the prior
(max-entropy when uniform) — flagged by a rank warning.

**Analytic λ SEs for `markov_ce` (Golan §4.7)** — `se_lambda` (S×K) uses the
**analytic Hessian** of the CE log-partition, `H_{(s,j),(s',j')} = Σ_k A_ks A_ks'
P_kj(δ_jj'−P_kj')` (eq. 4.7, the moment covariance under P̂; verified == the
finite-diff Hessian to machine precision, and now also the *bread* in
`.markov_sandwich` for the exact/no-noise case). **KEY subtlety:** `H⁻¹` *alone*
overstates the SE by **~√n_from** (measured ~22×) because `markov_ce`'s moments
`A`,`M` are cross-moment **sums** (`H = n·V̂`, so the information equality fails);
the reported SE is the robust **sandwich** `H⁻¹V̂H⁻¹` (Golan §3.3), which **matches
a unit bootstrap** (ratio ~1.0). λ is identified only up to a per-conditioning
additive shift (softmax), so `se_lambda` is **reference-normalized** to
destination `lambda_ref=1` (that column `NA`); a boundary fit (`P_kj∈{0,1}`, common
with covariates) → `NA` (infinite log-odds SE). Internals `.markov_hessian` +
`.markov_lambda_se`; object carries `se_lambda`/`lambda_ref`/`vcov_lambda`;
`summary()` prints it. See [[gme-se-validation]].

**Analytic λ SEs for `markov_gce` (Golan §4.7)** — the noisy sibling's dual has
**two** partition functions, so its analytic Hessian is the eq-4.7 covariance of
both: `(1/ν)·Σ_k A_ks A_ks' P_kj(δ_jj'−P_kj')` [signal] `+ δ_jj'·(1/(1−ν))·Σ_{i,t}
z1_its z1_its'·Var_w(v)_itj` [noise, **block-diagonal in destination j** because
each `Ψ_itj` depends only on column j of λ]. Unified in `.markov_info(sd, lhat)`
(reused as the exact **bread** in `.markov_sandwich`, replacing the finite-diff for
the noisy case — no accuracy change, just exact/faster; verified == finite-diff to
rel 5e-8). **KEY:** unlike `markov_ce`, the noise term **breaks the softmax
shift-redundancy**, so H is **full rank** and λ is fully identified → `se_lambda`
is a full S×K matrix (**no reference/NA column**). The reported SE is the
**unit-clustered** analytic sandwich `H⁻¹(Σ_unit g_unit g_unit')H⁻¹` (via
`.markov_lambda_se_full`), which **matches a unit bootstrap** (ratio ~1.0); the
non-clustered version overstates ~1.26× (covariate scores are correlated within a
unit). Object carries `se_lambda`/`vcov_lambda`; `summary()` prints it. **The
analytic Hessian does NOT improve the margins SEs** — the finite-diff bread it
replaces was already exact (5e-8), the margins AME already matches the bootstrap
(~1.1×), and clustering the *margins* meat makes it worse (1.57×, the AME
averaging washes out within-unit dependence), so margins is left non-clustered.
See [[gme-se-validation]].

**Inference add-ons (Golan §7.5/§7.7) for `markov_ce` + `markov_gce`** — the same
`fano_bounds()`/`margins()` generics as the multinomial/mixed estimators (shared
internals in `markov_ce.R`: `.fano_row_bounds` in `utils.R`, plus `.markov_margins`
/`.markov_sandwich`/`.markov_P_eps`/`.fd_jac`). (1) **`fano_bounds()`** — each row
of the K×K `p_hat` is a distribution over the next state, so this is **Golan's own
§7.7 Markov example**: per origin state `pe_k = 1 − max_j P_kj`, weak bound
`pe_k ≥ S(p_k) − ln2/ln(K)`; returns K rows (one per origin) + `"overall"`; a Fano
line is added to both `summary`/`print.summary`. Works for the plain no-covariate
`markov_ce` (the headline diagnostic there). (2) **`margins()`** — effect of a
covariate on the transition probabilities, `ME[k,j,s] = (P_kj/ν)(λ_sj − Σ_l P_kl
λ_sl)` (ν=1 for `markov_ce`), verified vs a numerical derivative to 7e-10.
**Requires covariates** (the no-covariate transition matrix has nothing to
differentiate → error). `average=TRUE` (default) → origin-frequency-weighted AME
(S×K); `average=FALSE` → the full K×K×S array. `se=TRUE` (with `average=TRUE`) attaches
SEs via `se_method = c("sandwich","bootstrap")` (default `"sandwich"`), returning
a `margins_gce` object. The naive Hessian delta overstates **~4.6×** vs a unit
bootstrap (median); the **sandwich** `(−H)⁻¹V̂(−H)⁻¹` (per-transition scores
`g_r,sj = z1_rs(y_rj−ε_rj) − z2_rs(x_r'P)_j`) **matches the bootstrap**
(~1.1–1.2× median), so it is the default. **`"bootstrap"`** resamples units at the
transition-block level using the stored `se_data` and **re-solves** the dual via
`.markov_refit` (standardization cancels, so it equals a refit-from-data
bootstrap — validated in tests); `se_data` carries a per-transition `unit` index.
Score inputs are stored in `object$se_data` (only when covariates present).
**Caveats:** (a) `markov_ce` with covariates is numerically finicky (the exact
covariate-conditioned moment problem often hits the optim cap — the reason
`markov_gce` exists); `margins` still computes on the returned P/λ, SE machinery
anchored on `markov_gce`. (b) For a **weakly-identified covariate** (e.g. a
near-deterministic column) the sandwich can *under*-state its SE where the
bootstrap is more honest — prefer `se_method="bootstrap"` there. See
[[gme-se-validation]].

**`markov_gce()`** — the noisy (GCE) sibling of `markov_ce()` (Golan §7.7.1, eq.
7.28). Same K×K row-stochastic Markov matrix, but each conditional moment carries
an additive error `ε_itj = Σ_m w_itjm v_m` on a symmetric `[−1,1]` support `v`,
weighted against the signal by `nu ∈ (0,1)`. This makes **over-identified**
moment systems — in particular **time-varying covariates** (which `markov_ce`'s
exact problem cannot bound) — feasible; **covariates are required**. **Reuses the
`.markov_arrays()` reshaper** from `markov_ce.R` (and therefore reads `a$Tn`, the
renamed list element). **By user instruction it keeps its own argument names**
(`data`, `id`, `time`, `states`, `covariates`, `v`, `nu`, `p0`, `w0`). **By user
choice the dual is kept in the maximise form (`fnscale=-1`)** — `markov_gce` and
`markov_ce` are the *two* maximise-form estimators (validated identical to the
minimise flip). Carries `p_hat` plus the noise as `epsilon` (Pn×K expected
per-observation errors) and `H_w` (no `w_hat` matrix); `H_signal = H_p` (mean row
entropy), prior-relative `S`. Two residual diagnostics: `foc_residual ≈ 0` (the
noisy first-order condition `M − A'P − Σ z ε`) vs `moment_residual` (the
exact-moment gap, *expected nonzero* — absorbed by the noise). As `v→0` or
`nu→0` it returns to `markov_ce`, but an **ultra-narrow `v` can diverge** (the
over-identified exact problem is unbounded — documented in `@section Choosing nu
and v`).

**`multinomial_gce()`** — the **nu-weighted GCE multinomial** of
Golan-Judge-Perloff (1996) with a **matrix `(y, X)` interface**. (Through v0.2.0
it had a formula-interface twin `gme_mnl()`, since removed; at `nu = 0.5`,
`which_alternative = 1`, and a matched support the two agreed exactly.) New parts:
the `nu` signal/noise weight, the matrix interface, and `which_alternative` (the
reference normalized to `λ = 0`). The uploaded draft had ~10 bugs (parse error,
undefined `e`/`w`, `v`/`mm` scoping, non-scalar objective, gradient dimension,
ones-vs-zeros reference, flat-par reshape) — all fixed; plus a data-scaled
`±1/√N` default support (the literal `±1` over-weights noise), log-sum-exp
stabilization, validation, and SEs of `lambda`. **By user instruction `coef()`
returns `lambda`** (the multipliers, K×J with the reference column 0), **not** a
`beta` (`beta = lambda/ν` would equal `gme_mnl`'s beta, but the user chose the
raw `lambda`). **Kept in the maximise form (`fnscale=-1`) by user choice** — the
third maximise-form estimator (with `markov_ce`/`markov_gce`). Adds a **new
exported `margins` generic** + `margins.multinomial_gce()` returning the K×J
average marginal effects `∂p_ij/∂x_ik = (p_ij/ν)(λ_kj − Σ_l p_il λ_kl)` (or the
N×J×K per-observation array via `average = FALSE`). Because the GCE noise forces
the reference alternative's error to zero, `which_alternative` slightly shifts
the fit (model property, not a bug). NOTE: the paper PDF is a scanned image
(no extractable text) — `gme_mnl` was the validation reference.

**Inference add-ons (Golan §7.5 / §3.6) for `multinomial_gce`** — two *distinct*
additions with *different provenance* (kept clearly separate in code + docs):
(1) **`fano_bounds()`** — a new exported generic + `fano_bounds.multinomial_gce()`
giving **Fano-inequality error bounds** for `p_hat` (an information-theoretic
*classification-error* bound, **not** a sampling SE). Per observation `i` (a row
of `p_hat`, a distribution over J): modal error `pe_i = 1 − max_j p_ij`, weak
bound `pe_i ≥ S(p_i) − ln2/ln(J)` with `S(p_i)=H(p_i)/ln(J)` (Golan's `−1/log K`
in bits = `−ln2/ln(J)` in nats). Returns a per-obs data frame (`p_max`/`pe`/`H`/
`S`/`pe_lower`) + an `"overall"` attribute; `summary()` prints a Fano line.
Validated: weak bound and the strong form `H(pe_i)+pe_i·ln(J−1) ≥ H(p_i)` both
hold for every obs. (2) **`margins(fit, se=TRUE)`** — SEs for the average marginal
effects, `se_method = c("sandwich","delta","bootstrap")`, **default `"sandwich"`**;
returns a shared `margins_gce` object with a `printCoefmat` table. **These SEs are
NOT from Fano** — Fano bounds classification error, these are sampling
variability. `"sandwich"` = the robust `(−H)⁻¹V̂(−H)⁻¹` with `V̂=Σ_i g_i g_iᵀ`
(per-obs free-λ scores = free-column block of `outer(X_i, Y_i−p_i−e_i)`)
sandwiched with the ME Jacobian — **matches the bootstrap** (ratio 0.83–0.96) and
is **PD/full rank** here (N≫K(J−1), no per-obs block, unlike `mixed_gce`).
`"delta"` = the classical naive `J·vcov·Jᵀ` (Golan §7.3 observed-information),
kept but **conservative ~1.35×** (the overstatement is `solve(−H)` under the GME's
finite-support regularization). Internal `.multinomial_gce_sandwich()`. NOTE —
not yet fixed: the object-level λ `se`/`vcov` (in `summary()`) are the same ~1.35×
naive Hessian; the sandwich is PD here so a `vcov_type` switch would fix them too
(easy follow-up). See [[gme-se-validation]].

**`mixed_gce()`** — the **doubly-reparameterized GCE "mixed" model** (from a
user math outline; Golan-style GCE). Response `Y` (N×J shares), design `X`
(N×J×K). Unlike `multinomial_gce` (which recovers `p` directly with a
reference-column normalization), here the signal probability is **itself**
reparameterized on a support `s ⊆ [0,1]`: `θ_ijm ∝ θ0_ijm exp(s_m(V_ij+ρ_i)/ν)`,
`V_ij = Σ_k X_ijk λ_kj`, `p_ij = Σ_m s_m θ_ijm`, with adding-up `Σ_j p_ij = 1`
enforced per-observation by a multiplier `ρ_i`; noise `w_ijh ∝ w0_ijh
exp(u_h V_ij/(1−ν))`, `e_ij = Σ_h u_h w_ijh`. Concentrated dual **maximized**
(`fnscale=-1`, with `markov_*`/`multinomial_gce`/`linreg_iv`/`panel_gce`) over
`(ρ, λ)`: FOCs are the adding-up `1−Σ_j p_ij=0` (in ρ) and the data moment
`Σ_i X_ijk(Y−p−e)=0` (in λ). **Genuinely new estimator.** **By user instruction
it keeps the draft's own argument names** (`Y`, `X`, `theta0`, `w0`, `nu`) and
**exposes `s`/`u`** (signal/noise supports; each accepts NULL / a point-count /
an explicit vector — `s ⊆ [0,1]`, `u` symmetric ⊆ `[−1,1]`). Uniform priors →
GME; user priors → GCE. **The uploaded draft never ran** — parse error in `grad2`
(two `colSums` with no operator, plus a dimensionally wrong `sweep(X,1,Y)`),
`est@par` (S4 on a list), the noise support `u` never passed to obj/grad, the
noise partition `Ψ` **missing `exp()` and the `w0` prior**, the objective ρ-term
`−ν·Σρ` (should be `+Σρ`; the `−ν` form diverges — `∂LL/∂ρ<0` always), and a
NULL-`w0` branch clobbering `theta0`. All fixed and validated: FOC/moment/
normalization ~1e-7, `p∈[0,1]`, valid `θ`/`w`, **shrink-to-pure recovery** of
`λ*/Pstar` as `u→0` (0.0002 at `u=±0.01`), GME=GCE-uniform, PD dual-Hessian SEs.
`coef()` returns `lambda` (K×J, GJP-family convention); carries `p_hat` and
`w_hat` plus `theta_hat` (N×J×M) and `rho_hat`/`se_rho`; `H_signal` is the
N-vector of per-obs signal (θ) entropies, `S = S_p`. **SEs** for λ and ρ from the
concentrated-dual Hessian (`vcov = solve(−H)`, validated PD). Adds
**`margins.mixed_gce`** on the existing `margins` generic (now dispatching for
both `multinomial_gce` and `mixed_gce`): `∂p_ij/∂x_ijk = (λ_kj/ν)·Var_θ(s)` —
the outline's eq. 8 **omitted the `1/ν`** (fixed; numeric check 3.5e-11).

**Inference add-ons (Golan §7.5 / §3.6) for `mixed_gce`** — mirrors the
`multinomial_gce` additions. (1) **`fano_bounds.mixed_gce`** — Fano error bounds
for `p_hat` (rows sum to 1 via ρ, so the row-distribution analysis applies
verbatim): `pe_i = 1 − max_j p_ij ≥ S(p_i) − ln2/ln(J)`; reuses the shared
internal `.fano_row_bounds(p, K)` (moved to `utils.R`; both estimators + both
summaries call it); `summary()` prints a Fano line. (2) **`margins(fit,
se=TRUE)`** — SEs for the average marginal effects. **KEY finding:** the naive
Hessian-inverse delta **overstates by ~7–9×** here (measured: λ's Hessian SE vs a
row bootstrap is 6.8–8.8×; the support reparameterization + per-obs ρ block make
`solve(−H)` a poor sampling-covariance estimate — far worse than
`multinomial_gce`'s 1.4×). So the default `se_method="sandwich"` uses the robust
`(−H)⁻¹V̂(−H)⁻¹` with `V̂=Σ_i g_i g_iᵀ` (per-obs dual scores: ρ-block `1−Σ_j p_ij`
≈0 at optimum → ρ pinned; λ-block `X_ijk(Y−p−e)`), which **matches the bootstrap
within ~13%** (validated); `se_method="bootstrap"` is the cross-check. The naive
Hessian delta is **not** offered for this model. Returns a shared `margins_gce`
object (the multinomial `margins_mnl_gce` class was renamed to this; single
`print.margins_gce`). **NOTE — not yet fixed:** `mixed_gce`'s object-level
`vcov`/`se_lambda`/`se_rho` (in `summary()`) are still the naive ~7–9× Hessian
inverse (the full sandwich is singular, rank ≤ N < N+K·J — a separate design
question); `summary()` now prints a caveat pointing to `margins(se=TRUE)` for
accurate ME SEs. See [[gme-se-validation]].

**`linreg_iv()`** — the **stochastic-moments GME estimator for instrumental-
variables regression** (Golan 2008, pp. 89–91), the IV sibling of `linreg`.
Model `y = Xβ + e` identified by instrument moments `IV'(y − Xβ − e) = 0`, with
`β = Zp` (signal support `Z`, K×M) and `e = vw`. Dual over the P instrument
multipliers `λ`: `p_km ∝ p0_km exp(z_km (X'IVλ)_k/ν)`, `w_ij ∝ w0_ij exp((IVλ)_i
v_j/(1−ν))`. **By user instruction it keeps its own argument names** (`y`, `X`,
`IV`, `Z`, `p0`, `v`, `w0`, `nu`). **The uploaded draft's math was wrong** — the
signal exponent used the *diagonal* `colSums(IV*X)` and the data term used `y'X`;
fixed to the **full `X'IVλ`** and **`y'IVλ`** (the gradient `IV'(y−Xβ−e)` was the
correct moment). Verified on the **David Card (1995)** return-to-schooling data
(n=3010): recovers no-intercept OLS (IV=X, `ed76 0.247`) and just-identified
2SLS (`ed76 0.139 ≈ 0.135`), where the draft gave `9.79` and diverged. Efficiency
improvements: precompute the constant `X'IV`/`IV'y`, log-sum-exp, and **internal
instrument standardization** (`cond(IV'X)≈2e6` → moment residual 414→2e-3, `β`
unchanged; `λ` reported on the original scale). **Generalised to over-identified
IV** (`λ` length = `ncol(IV) ≥ ncol(X)`). **Maximise form (`fnscale=-1`)** kept by
user choice (with `markov_*`/`multinomial_gce`). `coef()` returns `β`; carries
both `p_hat` and `w_hat`; `H_signal` is the K-vector of per-coef signal
entropies. **Standard errors are now implemented** via `se_method` ∈
{`sandwich` (default, robust HC0), `delta` (classical), `bootstrap` (pairs),
`none`}: `β = Zp` is a Z-estimator of the instrument moments, so
`Var(β) = Jβ A⁻¹ Ω A⁻¹ Jβ'` with `A` the dual Hessian,
`Jβ = ∂β/∂λ = diag(Var_p(z)/ν)·(X'IV)`, and `Ω` the meat (`IV'diag(r²)IV`
robust / `σ̂²·IV'IV` classical; `r = y−Xβ−e`). **MC-validated:** the sandwich
matches the trimmed MC sampling SD and **→ 2SLS robust SE as the support widens**
(0.1069 vs 0.1068), the classical delta assumes homoskedasticity (understates
under heteroskedasticity). Adds `vcov.linreg_iv(type=)` + an ER-style
`summary()` coef table; object carries `se_beta`/`vcov`/`se_method`/`boot`. SEs
are **asymptotic** (the β sampling distribution is support-bounded/skewed — Wald
intervals approximate). See [[gme-se-validation]]. NOTE: Card's paper PDF is a
scanned image; OLS/2SLS from the data were the validation reference.

**`inverse_noise()`** — the noisy-moment (GME/GCE) sibling of `inverse_ce()`,
solving `y = X p + e` with `e_t = Σ_j v_j w_tj` (Golan §6.1). Same lm-style
formula orientation (response = moments, RHS = states). **By deliberate user
instruction it keeps its own argument names** (`p0` signal prior, `w0` T×M noise
prior, `v` noise support, `nu` ∈ (0,1) signal/noise entropy weight) rather than
`me()`'s `q` — and these are in fact the standard GME/GCE names from this file's
conventions. Uniform `p0`/`w0` → GME; user priors → GCE. It minimises the convex
dual `-λ'y + (1−ν)log Ω + ν Σ log Ψ_t` (no `fnscale`), matching `gme()` and
`inverse_ce`. Carries both `p_hat` and `w_hat`. As `v → 0` the estimates converge
to the `inverse_ce` solution.

**De-duplicated output + dropped `converged` (user request).** The object keeps
only the canonical member of each former duplicate pair: `p_hat`, `lambda_hat`,
`w_hat`, **`residuals`** (kept over `e` by user choice — the moment residual
`y−Xp` *is* the estimated noise `e` at the optimum), and `objective`; the aliases
`p`/`lambda`/`w`/`e`/`value` are gone. **`converged` is dropped, only the raw
`convergence` code is kept** — a deliberate **exception** to the "every object
carries `converged`" canonical-minimum convention in this file, for
`inverse_noise` only. (Aside: the noise-prior entropy field is named `H_error0`,
not `entropy_w0`, to avoid `$e` partial-matching it now that `e` is gone.) The
solver is factored into **`.inverse_noise_engine(X, y, v, nu, p0, w0, con)`**
(reused by the bootstrap). Sibling `inverse_ce` still carries its duplicates — not
touched (out of scope).

**Sampling SEs for `lambda_hat` and `p_hat` (validated by MC).** Structurally
`inverse_noise` is a GME regression `y=Xp+e` with `p` on the simplex, so the **T
moment rows are the sample size** and sampling SEs are meaningful. Since `λ̂`
solves `Xp(λ)+e(λ)−y=0`, the IFT gives `∂λ̂/∂y=H⁻¹`, so `Var(λ̂)=H⁻¹Σ_e H⁻¹` and
`Var(p̂)=J·Var(λ̂)·Jᵀ` with `J_{k,t}=∂p_k/∂λ_t=(1/(1−ν))p_k(x_tk−(Xp)_t)`. Three
methods via `se_method=` (fit-time, default **`"sandwich"`**), `vcov(type=)`, and
`summary(se_method=)`: **`sandwich`** = robust HC0 `H⁻¹ diag(residuals²) H⁻¹`
(meat = squared fitted noise) — the accurate default (validated ≈ MC SD and
bootstrap; per-element `λ` SEs are noisier, each resting on one `ê_t²`);
**`delta`** = naive `H⁻¹` — kept for comparison but **overstates ~8–11×** (NOT a
sampling covariance); **`bootstrap`** = residual bootstrap (rows fixed, resample
`ê`, refit via the engine, `B` reps) — stable aligned SEs for **both** `λ` and `p`,
the gold cross-check. **`vcov()`'s default changed from the raw `H⁻¹` to the
sandwich** (`type="delta"` recovers the old `H⁻¹`); all three mildly **understate
(~15–20%)** because the finite-support entropy penalty shrinks `ê` below the true
noise dispersion (the bootstrap shows the same — estimator property, not a bug).
Object stores `se_lambda`/`se_p`/`vcov_lambda`/`se_method` (+ `X`/`y`/`control` for
recompute/bootstrap; `NULL` SEs when `se_method="none"`). See [[gme-se-validation]].

**Signal + noise normalized entropy and Fano (user request).** Reports **both**
the signal `S = S_p = H(p̂)/H(p0)` (existing) and the **noise** `S_w = H(ŵ)/H(w0)`
(new, mirroring `matrix_gce`'s `ent(M)=−ΣM log M`; fields `H_error`/`H_error0`);
`print`/`summary` show `signal S(p) …, noise S(w) …`. **`fano_bounds.inverse_noise`**
is the verbatim parallel of `fano_bounds.inverse_ce` (the recovered `p̂` is one
distribution over K states): `.fano_row_bounds(matrix(p̂,1,K),K)`, `pe=1−max_k p̂_k
≥ S(p̂)−ln2/ln K`; `summary()` prints a Fano line. `summary()`'s `p` table is
`Estimate | Std. Error | t value` (via `printCoefmat`).

**`matrix_gce()`** — the noisy (GCE) sibling of `matrix_ce()` (Golan §7.4),
solving `y = Px + e` with `e_i = Σ_h v_h w_ih`. Recovers a **column-stochastic**
signal `p` (n×m, columns sum to 1) **and** a **row-stochastic** noise `w` (n×H,
rows sum to 1) via the dual over the n row multipliers:
`p̂_ij ∝ p0_ij exp(λ_i x_j/ν)`, `ŵ_ih ∝ w0_ih exp(λ_i v_h/(1−ν))`. **By user
instruction it keeps its own argument names** (`y`, `x`, `p0` signal prior,
`w0` noise prior, `v` noise support, `nu` ∈ (0,1) entropy weight). `y`/`x` are
normalised to `[0,1]` so errors lie in `[−1,1]` (warns otherwise); `v` defaults
to `c(−1,0,1)`. Uniform priors → GME; user priors → GCE. Minimises the convex
dual `−Σ λ_i y_i + ν Σ_j log Ω_j + (1−ν) Σ_i log Ψ_i` (no `fnscale`, flip
validated against the maximise form — identical `p` and `w`), matching `gme()`/
`inverse_noise`. Carries both `p_hat` and `w_hat`; `H_signal` is the scalar
`H(p̂)`, the canonical `S` is the signal normalized entropy `S_p = H(p̂)/H(p0)`
(with `S_w` alongside). `summary()` reports the Golan §7.5 information measures
(normalized entropy `S(P)`, information index `I = 1−S`, pseudo-R², per-column
`S(p_j)`) and an entropy-ratio test `W = 2[H(p0)−H(p̂)] ~ χ²(n−1)` for
H0: P = P0 — **preserved verbatim from the uploaded code**. As `v → 0` the
signal converges to the `matrix_ce` solution.

**`linreg()`** — fits `y = Xβ + e` by GCE (the concentrated dual, `nu` weight,
and Golan p. 96 asymptotic `vcov`) with an lm-style interface: a `predict()`
method, `r.squared`, and an **entropy-ratio (ER) `summary()`**. (Through v0.2.0
this shared its estimator with the now-removed `gme()`.) The ER summary reports, per coefficient, an ER test of
H0: β_k = 0 (Golan §6.4/§6.6): `ER_k = 2[H*_unrestricted − H*_restricted,k]`
with `H* = Σ H(p) + Σ H(w)`, restricted fit zeroing row k of the signal support
`Z`, `ER_k ~ χ²(1)`; plus an overall test of H0: all coefficients (intercept and
slopes) = 0 (zero **every** `Z` row; `ER ~ χ²(K)`). The coefficient table is
`Estimate | Std. Error | ER | Pr(>Chi)` (via `printCoefmat`). Mechanics: the
inline solver was factored into the internal engine `.linreg_engine(X, y, Z, p0,
v, w0, nu, con, init)` (the main fit and the K+1 warm-started restricted refits
share it) plus `.linreg_Hstar(p, w)`; ER can be slightly negative (collapsing a
support row raises that coef's entropy to its prior max) so it is clamped at 0.
To make refits faithful, `linreg()` **stores the resolved `Z`, `p0`, `v`,
`w0`, `X`, `y`, `control`** on the object. Its **argument names** use lowercase
`v` for the noise support, and otherwise follow the package conventions (`Z`,
`p0`, `w0`, `nu`; uniform priors → GME, user priors → GCE; minimised convex
dual; PD Hessian).
`linreg()` carries both `p_hat` (K×M) and `w_hat` (T×J); its `H_signal` is the
K-vector of per-coefficient signal entropies, and its normalized entropy `S` is
**prior-relative** — `S = H(p̂)/H(p0)` (Golan §6.4), reducing to `H/(K·log M)`
for a uniform prior (GME) and matching `matrix_ce`/`markov_ce`/`inverse_noise`;
it can exceed 1 for a strongly informative prior.

**`panel_gce()`** — dual **GME/GCE for the one-way error-components panel model**
`y_nt = x_nt'β + μ_n + ε_nt` (Lee & Cheon 2014, *CSAM* 21(5), eq. 3.12), the
panel/individual-effects sibling of `linreg`. Three reparameterized
blocks are recovered jointly: the coefficients `β_k = Σ_m z_km p_km` (support
`Z` K×M, prior `p0`), the **individual effects** `μ_n = Σ_r f_nr g_nr` (support
`FF` N×R, prior `g0`), and the errors `e_nt = Σ_j v_j w_ntj` (support `v`, prior
`w0`), via the concentrated dual over `λ` (length N·T): `p_km ∝ p0_km exp(z_km
(X'λ)_k/ν)`, `g_nr ∝ g0_nr exp(f_nr (Σ_t λ_nt)/ν)`, `w_ntj ∝ w0_ntj exp(v_j
λ_nt/(1−ν))`; the FOC is the model equation `y − Xβ − μ − e = 0`. The `nu`
weights signal (β, μ) vs noise; default `v` is a 3σ grid on the within-residual
sd (Swamy-Arora spirit). **Two `layout`s** map `λ` to `(n,t)`: `period` (index
`(t−1)N+n`, units stacked within period) vs `unit` (index `(n−1)T+t`, periods
stacked within unit) — set to match how the panel rows are ordered.
**The draft's math was verified correct against eq. 3.12 (no bugs)**; on a
simulated balanced panel it recovers the within/fixed-effects estimate
(`max|β̂ − β_FE| ≈ 0.005`, both layouts) and the individual effects
(`cor(μ̂,μ) ≈ 0.98`). **By user instruction it keeps its own argument names**
(`y`, `X`, `Z`, `tt`, `FF`, `p0`, `g0`, `v`, `w0`, `nu`, `layout`). **Maximise
form (`fnscale=-1`)** kept by user choice — the paper minimises the `λ→−λ`
equivalent (with `markov_*`/`multinomial_gce`/`linreg_iv`). Uniform `p0`/`g0`/`w0`
→ GME; user priors → GCE. **As in `markov_ce`, β is prior-insensitive when
well-determined** (a wide `Z` makes the prior nearly irrelevant; the prior bites
only on a moderate/narrow support — the test checks the *directional* pull, not a
magnitude threshold). Carries `p_hat`/`g_hat`/`w_hat`, `mu_hat`, `e_hat`,
`H_p`/`H_g`/`H_w`; `H_signal` is the K-vector of per-coefficient signal
entropies, `S = S_p` the normalized signal entropy. **Standard errors (Lee &
Cheon §3.3) are now implemented** via a `vcov_type` argument, MC-validated
(β 0.7%, p 1%, μ 1.2%): **β** — `"within"` (default) `σ̂²_ε(X̃'X̃)⁻¹`, the
*accurate* SE because `μ̂` makes β within-identified (the paper's raw-`X`
sandwich overstates by ~15–20%); `"cluster"` (within, unit-clustered robust);
`"disturbance"` (the literal §3.17 form `(X'X)⁻¹ X'Σ̂_uX (X'X)⁻¹`, conservative).
**p** — delta from β: `SE(p̂_km)=|p_km(z_km−β_k)|/Var_p(z_k)·SE(β̂_k)`. **μ** —
finite-T prediction SE `√(σ̂²_ε/T_n + x̄_n'Var(β̂)x̄_n)` (individual effects are
not √N-consistent). Object carries `vcov` (K×K, via `vcov()`), `se_beta`,
`se_p` (K×M), `se_mu` (N), `sigma2_eps`, `vcov_type`; `summary()` prints a coef
table with z & Pr(>|z|). All three SE formulas use only X/residuals/`p_hat`/`β̂`,
so they are **independent of the `fnscale=-1` sign convention**. **X must not
include an intercept** (μ subsumes the level; a constant column makes the within
design singular). Note: the dual can rarely (~0.3%) converge to a spurious
optimum from the `λ=0` start (bad basin) — the analytic SEs describe the good
optimum.

Every estimator object must contain at minimum:
`p_hat`, `lambda_hat`, `H_signal`, `S`, `objective`, `converged`, `method`, `call`.

`inverse_ce()`, `matrix_ce()`, and `markov_ce()` carry
`p_hat` (no `w_hat`). `inverse_noise()`, `linreg()`, `matrix_gce()`,
`multinomial_gce()`, `linreg_iv()`, `panel_gce()`, and `mixed_gce()` carry both
`p_hat` and `w_hat` (for `multinomial_gce`, `coef()` returns `lambda`, not a
`beta` — by user instruction; for `linreg_iv`, `coef()` returns `β = Zp`, and
`lambda_hat` is the P-vector of instrument multipliers; for `panel_gce`, `coef()`
returns `β = Zp`, and it additionally carries `g_hat` — the individual-effects
weights, with `mu_hat = F·g` — and `lambda_hat` is the N·T-vector of dual
multipliers; for `mixed_gce`, `coef()` returns `lambda` (K×J) and it additionally
carries `theta_hat` (N×J×M, the doubly-reparameterized signal weights with
`p_hat = Σ_m s_m θ`) and `rho_hat`/`se_rho` — the per-observation adding-up
multipliers). `markov_gce()`
carries `p_hat` plus the noise as `epsilon` (Pn×K expected per-observation
errors) and `H_w` — no `w_hat` matrix (the noise weights form a Pn×K×M array;
`epsilon` is the 2D summary). For
`inverse_ce()` (de-duplicated, same cleanup as `inverse_noise`), `H_signal` is the
scalar `H(p_hat)` and `objective` is the minimised dual value; it carries
`p_hat`/`lambda_hat`, `fitted.values`, `residuals`, `p0` (the ME/CE prior — the
arg is `p0`, renamed from `q`), `hessian`, and the information-matrix inference
fields `vcov_lambda`/`se_lambda`/`se_p` (`NA`/`NULL` when `I(λ)` is singular); it
**drops** `converged` (keeps `convergence`) and the `p`/`lambda`/`value`/`prior`/`q`
aliases. `inverse_noise()`
(de-duplicated) carries `p_hat`/`lambda_hat`/`w_hat`, `fitted.values`,
`residuals` (= the noise `e`), `objective`, `nu`, `prior`/`noise_prior`,
`support` (the `v` vector), the entropy fields `H_signal`/`H_error`/`H_error0`
and `S`/`S_p`/`S_w`, the SE fields `se_lambda`/`se_p`/`vcov_lambda`/`se_method`,
and `X`/`y`/`control` (for SE recompute/bootstrap); it **drops** `converged`
(keeps `convergence`) and the `p`/`lambda`/`w`/`e`/`value` aliases.

### Support spaces
Support spaces Z (for β) and V (for ε) are the user's most consequential
input in GME/GCE. Always:
1. Validate that the support is symmetric around zero for V (error support).
2. Validate that the support for each β_k contains at least M=2 points.
3. The three-sigma rule (`error_scale = 3`) is the default for V.
4. Document clearly that wider supports = more regularization toward uniform.

---

## Code conventions

### Naming
- **Functions**: `snake_case`, descriptive verbs for estimators (`me`, `gme`,
  `el`, `gel`), nouns for measures (`shannon_entropy`, `kl_divergence`).
- **Internal helpers**: prefix with `.` (e.g., `.me_dual`, `.check_prob`).
  These are not exported and need no `@export` tag.
- **Arguments**: `y` for the dependent/moment vector, `X` for the design/moment
  matrix, `q` for ME/CE prior, `p0` for GME/GCE signal prior (K × M matrix),
  `w0` for GME/GCE error prior (T × J matrix), `Z` for signal support,
  `V` for error support, `lambda` for Lagrange multipliers, `p` for
  probabilities.
- **No abbreviations** in public argument names (`half_range` not `hr`,
  `method` not `meth`).

### Documentation
Every exported function must have:
- `@param` for every argument
- `@return` describing the output precisely
- `@details` with the mathematical formula and its equation number from
  Golan (2008) where applicable
- `@references` citing at minimum Golan (2008) and the original source
- `@examples` that run in under 5 seconds with no external data dependencies
- `@export`

Internal helpers need only a one-line `#' @keywords internal` comment.

### Error messages
Follow this pattern — always state what was received:
```r
stop("X must have nrow(X) == length(y). Got nrow(X) = ", nrow(X),
     ", length(y) = ", length(y), ".")
```
Never use `T` or `F` as abbreviations for `TRUE` / `FALSE` (CRAN policy).

### No global state
No `<<-`, no `options()` side-effects, no `library()` or `require()` calls
inside any function. Use `requireNamespace("pkg", quietly = TRUE)` to test
for optional packages.

---

## Testing conventions

Tests live in `tests/testthat/`, one file per source file:
`test-entropy.R`, `test-utils.R`, `test-linreg.R`, `test-inverse_ce.R`, etc.

### What every test file must cover
1. **Known analytical results** — test against closed-form truths (e.g.,
   `shannon_entropy(uniform) == log(K)`, `kl_divergence(p, p) == 0`).
2. **Limiting cases** — verify that special cases reduce correctly
   (e.g., CE with uniform prior == ME, EL with OLS moments == OLS,
   **GCE with uniform priors == GME**).
3. **Moment constraint satisfaction** — `max(abs(residuals(fit))) < 1e-4`
   at convergence for dual solver.
4. **Input validation** — confirm that malformed inputs throw the right errors.
5. **S3 methods** — confirm `print`, `summary`, `coef`, `fitted`, `residuals`
   run without error and return objects of the right length.

### Tolerance policy
The dual BFGS solver achieves moment residuals of roughly 1e-5 to 1e-6 for
well-scaled problems. Use `tolerance = 1e-4` for moment satisfaction checks,
`tolerance = 1e-8` for exact algebraic identities (e.g., CR reduces to KL).

---

## Key mathematical reference (Golan 2008)

| Concept | Equation | Location |
|---|---|---|
| Shannon entropy | H(p) = −Σ p_k log p_k | Eq. 3.1 |
| KL divergence | D(p‖q) = Σ p_k log(p_k/q_k) | Eq. 3.4 |
| Rényi entropy (order α) | H^R_α = (1/(1−α)) log Σ p^α_k | Eq. 3.5 |
| Tsallis divergence (order α) | D^T_α = (1/(1−α))(Σ p^α q^{1−α} − 1) | Eq. 3.7 |
| Cressie-Read (order α) | D^CR_α = (1/(α(1+α))) Σ p((p/q)^α − 1) | Eq. 3.8 |
| Unifying CR relationship | D^R_{α+1} = −(1/α)log[1+α(α+1)D^CR_α] | Eq. 3.9 |
| ME primal | max H(p) s.t. Xp = y, Σp = 1 | Eq. 4.1 |
| CE primal | min D(p‖q) s.t. Xp = y, Σp = 1 | Eq. 4.2 |
| CE solution | p̃_k = q_k exp(λ'x_k) / Ω(λ) | Eq. 4.3 |
| Dual concentrated | min_λ {−λ'y + log Ω(λ)} | Eq. 4.4–4.5 |
| Normalized entropy | S(p̃) = H(p̃) / H(p⁰) | Section 6.4 |
| Entropy-ratio statistic | W = 2K log(K)(1 − S) ~ χ²(K−1) | Section 6.4 |
| Pseudo-R² | 1 − S(p̃) | Section 7.5 |
| GME primal | max H(p)+H(w) s.t. y=XE_P[Z]+E_W[V] | Eq. 6.5 |
| GCE primal | min D(p‖p0)+D(w‖w0) s.t. y=XE_P[Z]+E_W[V] | Eq. 6.5 (general) |
| GME/GCE dual | min_λ {Σy_t λ_t + Σ log Ω_k(λ) + Σ log Ψ_t(λ)} | Below Eq. 6.9 |
| GCE signal probs | p̂_km = p0_km exp(−z_km Σ λ_t x_tk) / Ω_k(λ̂) | Eq. 6.6 (general) |
| GCE error probs | ŵ_tj = w0_tj exp(−λ̂_t v_j) / Ψ_t(λ̂) | Eq. 6.7 (general) |
| GME signal probs | p̂_km = exp(−z_km Σ λ_t x_tk) / Ω_k(λ̂) | Eq. 6.6 (p0 uniform) |
| GME error probs | ŵ_tj = exp(−λ̂_t v_j) / Ψ_t(λ̂) | Eq. 6.7 (w0 uniform) |

---

## CRAN compliance checklist (enforce before every commit)

- [ ] No `T` / `F` as logical literals
- [ ] No `library()` / `require()` inside functions
- [ ] No `<<-` global assignment
- [ ] All exported functions have `@examples` running in < 5 seconds
- [ ] `NAMESPACE` matches actual exports (regenerate with `devtools::document()`)
- [ ] `R CMD check` produces zero ERRORs and zero WARNings
- [ ] `Description:` field in `DESCRIPTION` is two or more complete sentences
  and contains no URLs
- [ ] Package title does not start with "A" or "The"
- [ ] For `License: GPL-3`, do NOT ship a `LICENSE` file (CRAN provides the
  standard GPL text; a bundled copy triggers a NOTE)
- [ ] `BugReports:` and `URL:` fields populated in `DESCRIPTION`

---

## Do not

- Do not add `Rcpp` or C/C++ code without discussion — the package is
  designed to be pure R for simplicity and portability.
- Do not change the default solver from BFGS dual to anything else without
  benchmarking and updating this file.
- Do not hard-code support spaces inside estimator functions — always expose
  them as user arguments with documented defaults.
- Do not use `set.seed()` inside any exported function.
- Do not write vignette code that requires internet access or takes more than
  30 seconds to run.
- Do not create separate `gce_*` twins of the GME estimators. Each estimator
  subsumes GCE via its `p0`/`w0` prior arguments (uniform priors → GME). One
  function per estimator family.
