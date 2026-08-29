# CRAN submission comments — infometrics 0.3.0

This is a new submission.

## Test environments

* Local: Windows 11 x86_64-w64-mingw32, R 4.5.2
* win-builder: R-devel  — **TODO: run and record the result before submitting**
* win-builder: R-release — **TODO: run and record the result before submitting**

## R CMD check results

Local (`R CMD check --as-cran`, R 4.5.2, Windows 11):

```
0 errors | 0 warnings | 0 notes
```

On CRAN's incoming checks this being a first submission is expected to produce
the usual NOTE:

```
Maintainer: 'Ganbaatar Jambal <jg3169a@gmail.com>'
New submission
```

## Notes for the reviewer

* The package is pure R with no compiled code, and imports only `stats`.
* The reference in the `Description` field is written in the required form:
  Golan (2008) <doi:10.1561/0800000004>.
* All exported functions are documented with `\value` sections, and all
  examples run unwrapped (no `\dontrun{}` or `\donttest{}`).
* The package writes nothing to the user's home directory, the working
  directory, or the package library; no examples, tests, or vignettes access
  the internet.
* Two vignettes are included and are rebuilt at check time (both use only
  simulated data and run in a few seconds each).
* This is a new package with no reverse dependencies, so no reverse-dependency
  checks were required.

## Test suite

689 `testthat` tests pass with 0 failures and 0 warnings.
