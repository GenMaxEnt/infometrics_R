# infometrics

**Information-Theoretic Methods for Econometric Estimation**

`infometrics` is an R package implementing the class of Information-Theoretic
(IT) estimators for econometric models, following the unified framework of
Golan (2008), *Information and Entropy Econometrics — A Review and Synthesis*
(Foundations and Trends in Econometrics, 2(1–2), 1–145). The central unifying
concept is the Cressie–Read power-divergence criterion, a generalized entropy of
order α that nests Shannon entropy, Empirical Likelihood, and the Pearson
chi-squared statistic as special cases.

All estimators solve the **dual (concentrated) model** by default — working in
Lagrange-multiplier space via BFGS with analytic gradients — which is
computationally efficient and needs no compiled code (the package is pure R,
depending only on `stats`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("GenMaxEnt/infometrics_R")
```

(The repository is currently private; installing requires access and a
GitHub token, e.g. `remotes::install_github("GenMaxEnt/infometrics_R", auth_token = ...)`.)

## What's included

**Entropy & divergence measures** — Shannon entropy, Kullback–Leibler divergence,
Rényi entropy/divergence, Tsallis divergence, Cressie–Read, normalized entropy.

**Maximum Entropy / Cross-Entropy (ME/CE)**

- `me()` — ME/CE estimator (matrix/vector interface), primal and dual.
- `inverse_ce()` — pure-moment ME/CE inverse problem (formula interface), with
  information-matrix standard errors and Fano error bounds.
- `matrix_ce()` — CE/ME matrix balancing.
- `markov_ce()` — CE Markov transition-matrix estimation.

**Generalized ME / Generalized CE (GME/GCE)** — one function per family;
uniform priors give GME, user priors give GCE.

- `gme()` — GME/GCE linear regression; `linreg()` — lm-style GCE regression with
  an entropy-ratio summary; `linreg_iv()` — GME instrumental-variables regression.
- `inverse_noise()` — noisy-moment GME/GCE inverse problem, with sandwich /
  delta / bootstrap standard errors and Fano bounds.
- `matrix_gce()` — noisy matrix balancing; `markov_gce()` — noisy Markov chains.
- `me_mnl()`, `gme_mnl()`, `gce_clogit()`, `multinomial_gce()`, `mixed_gce()` —
  multinomial and (mixed / conditional) logit estimators.
- `panel_gce()` — one-way error-components panel model.

Every estimator returns an S3 object inheriting from `"infometrics"` with
`print`, `summary`, `coef`, `fitted`, and `residuals` methods (and, where
applicable, `vcov`, `margins`, and `fano_bounds`).

## Reference

Golan, A. (2008). *Information and Entropy Econometrics — A Review and
Synthesis*. Foundations and Trends in Econometrics, 2(1–2), 1–145.

## License

GPL-3.
