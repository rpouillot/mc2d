## Submission

This is an update of mc2d (0.3.0), following a full review of the code:

* Bug fixes (see inst/NEWS).
* Faster code and lower memory use on large mcnodes (Monte-Carlo arrays).
  The results obtained with a given seed are unchanged, except for the
  corrected bugs (listed in inst/NEWS).
* The dependency on 'ggpubr' (and its ~60 recursive dependencies) was removed:
  the arrangement of several ggplot2 graphs is now done internally with
  'gtable' (already a dependency of 'ggplot2') and 'grid'.
* No change in the API: no exported function, argument or default value was
  removed or modified.
* Many new tests (tests/testthat).

## Test environments

* local Windows 11, R 4.6.0 (R CMD check --as-cran): Status OK
* local Windows 11, R-devel (2026-10-20 r90728 ucrt) (R CMD check --as-cran): Status OK
* win-builder (R-devel): Status OK
* R-hub (R-devel: linux, windows, macos, macos-arm64): Status OK

## R CMD check results

0 errors | 0 warnings | 0 notes

(local, win-builder and R-hub checks)

## Reverse dependencies

We checked the 20 reverse dependencies (aae.pop, bayespmtools, CloneSeeker,
CopulaREMADA, CPP, decisionSupport, fitdistrplus, fitteR, FitVerse, handwriter,
hmix, mcmodule, OSCV, PRA, rriskDistributions, sensitivity, stratifyR, Umpire,
wired, wxgenR) with R CMD check (local Windows 11, R 4.6.0, 2026-10-22),
comparing the current CRAN version of mc2d (0.2.2) with this version (0.3.0).

 * We saw 0 new problems: the check results are identical with both versions.
 * We failed to check 0 packages.
