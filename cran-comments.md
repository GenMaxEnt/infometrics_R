# CRAN submission comments — infometrics 0.3.0

This is a new submission.

## Test environments

* Local: Windows 11 x86_64-w64-mingw32, R 4.5.2
* win-builder, x86_64-w64-mingw32: R-devel (2026-08-27 r90452)
* win-builder, x86_64-w64-mingw32: R release 4.6.1 (2026-06-24)

## R CMD check results

Local (`R CMD check --as-cran`, R 4.5.2, Windows 11):

```
0 errors | 0 warnings | 0 notes
```

Both win-builder runs (R-devel and R-release) returned `Status: 1 NOTE`, the
same NOTE in each case:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Ganbaatar Jambal <jg3169a@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  Fano (18:24)
  GCE (13:6)
  GME (12:43)
```

Both parts of this NOTE are expected and, I believe, benign:

* **New submission** — this is the first submission of the package to CRAN.
* **Possibly misspelled words** — all three are spelled correctly and are
  standard terminology in this literature:
  * **GME** and **GCE** are the established acronyms for Generalized Maximum
    Entropy and Generalized Cross-Entropy. Both are written out in full at
    first use in the `Description` field, immediately before the acronym.
  * **Fano** is a surname (Robert M. Fano), used here for Fano's inequality;
    the package reports Fano error bounds through `fano_bounds()`.

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
