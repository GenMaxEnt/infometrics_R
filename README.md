# infometrics

**Information-Theoretic Methods for Econometric Estimation**

`infometrics` is an R package implementing Information-Theoretic (IT) estimators
for econometric models, following the unified framework of Golan (2008),
*Information and Entropy Econometrics — A Review and Synthesis*
([doi:10.1561/0800000004](https://doi.org/10.1561/0800000004)). The package
focuses on the Generalized Maximum Entropy (GME) and Generalized Cross-Entropy
(GCE) family: each unknown is reparameterized as the mean of a distribution on a
bounded support, and the resulting entropy problem is solved subject to the data.

Every estimator solves the **dual (concentrated) model** — working in
Lagrange-multiplier space via BFGS with analytic gradients — which is
computationally efficient and needs no compiled code. The package is pure R and
depends only on `stats`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("GenMaxEnt/infometrics_R")
```

## What's included

One function per estimator family: **uniform priors give GME, user-supplied
priors (`p0`, `w0`) give GCE** — there are no separate `gce_*` twins.

**Regression**

- `linreg()` — lm-style GCE linear regression, with `predict()`, an
  entropy-ratio (ER) coefficient test, and prior-relative normalized entropy.
- `linreg_iv()` — GME instrumental-variables regression (just- and
  over-identified), with sandwich / delta / bootstrap standard errors.
- `panel_gce()` — one-way error-components panel model, with standard errors
  for the coefficients, individual effects, and support weights.

**Inverse problems**

- `inverse_ce()` — pure-moment ME/CE inverse problem (formula interface), with
  information-matrix standard errors and Fano error bounds.
- `inverse_noise()` — noisy-moment GME/GCE inverse problem, with sandwich /
  delta / bootstrap standard errors and Fano bounds.

**Matrix balancing and Markov chains**

- `matrix_ce()` / `matrix_gce()` — exact and noisy-moment matrix balancing.
- `markov_ce()` / `markov_gce()` — exact and noisy-moment first-order Markov
  transition matrices, with marginal effects and Fano bounds.

**Multinomial response**

- `multinomial_gce()` — nu-weighted GCE multinomial (Golan–Judge–Perloff),
  with marginal effects (robust SEs) and Fano bounds.
- `mixed_gce()` — doubly-reparameterized GCE mixed model.

**Utilities** — `shannon_entropy()`, `make_support()`, `default_supports()`,
`normalize_data()`, plus the `margins()` and `fano_bounds()` generics.

Every estimator returns an S3 object inheriting from `"infometrics"` with
`print`, `summary`, `coef`, `fitted`, and `residuals` methods (and, where
applicable, `vcov`, `predict`, `margins`, and `fano_bounds`).

## Vignettes

```r
vignette("inverse-problems",    package = "infometrics")
vignette("gme-gce-regression",  package = "infometrics")
```

- *Pure Inverse Problems* — maximum-entropy and cross-entropy recovery of a
  distribution from moment constraints.
- *Generalized Maximum Entropy and Cross-Entropy for the Linear Model* —
  `linreg()` and `linreg_iv()`: GME under collinearity, GCE with an informative
  prior, the entropy weight `nu`, entropy diagnostics and hypothesis tests, and
  an instrumental-variables example.

## A note on support spaces

Support spaces are the user's most consequential input. The dual is bounded
below only if the data are representable within the supports — if they are not,
the multipliers diverge and estimates collapse to a support boundary. The
estimators check the first-order condition at the optimum (`foc_residual`) and
warn when it is violated; the remedy is a wider error support `v` (and/or `Z`).

## Reference

Golan, A. (2008). *Information and Entropy Econometrics — A Review and
Synthesis*. Foundations and Trends in Econometrics, 2(1–2), 1–145.
[doi:10.1561/0800000004](https://doi.org/10.1561/0800000004)

## License

GPL-3.
