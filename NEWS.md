# cit 3.0.0

## New features

* `cit.effects()` estimates the total, direct and indirect effects of L on T through G,
  and the proportion mediated, with delta-method confidence intervals, for continuous
  (mean difference), binary (odds ratio) and survival (hazard ratio) outcomes. It uses
  the published regression estimators without an exposure-mediator interaction
  (VanderWeele and Vansteelandt 2009, 2010; VanderWeele 2011; Valeri and VanderWeele
  2013); the ratio scales assume a rare outcome. Results agree with the `regmedint`
  package. The estimates need stronger assumptions than the CIT (no unmeasured
  confounding of G and T); see `?cit.effects`.

* Sign-pattern CIT (`sign.pattern = TRUE`, the default in `cit.cp()` and `cit.bp()`
  with a single instrument and a single mediator). Mediation now also requires the
  estimates in tests 1-3 to agree in sign, sign(T ~ L) = sign(T ~ G | L) x
  sign(G ~ L | T). The alternative is the union of the four sign patterns that
  satisfy this; within the pattern supported by the data, tests 1-3 are one-sided.
  Associations with the wrong sign (for example from a hidden confounder of G and T
  that reverses test 3) no longer count as evidence for mediation, and when the
  signs agree, tests 1-3 report half their two-sided p-values. The error rate is
  controlled over the whole sign-null (proof for the test 2 and test 3 null faces,
  simulation elsewhere). With `robust = TRUE` the two-sided p-values are kept and a
  component with the wrong sign is set to 1. `sign.pattern = FALSE` gives the
  two-sided components of 2.3.2.

* `cit.cox()` for a censored survival outcome. Test 3 is the test of G given the
  event indicator (L in G ~ E + CG + L), which guards against the event changing G
  (for example, tissue collected at or after the event). `test3 = FALSE` omits it
  when G was measured before any event could occur (see `?cit.cox`).
* Multiple mediators: `G` may be a matrix of mediators tested jointly (test 3 uses
  Wilks' lambda).
* Separate covariates for the outcome models (`CT`, formerly `C`) and for the
  mediator model (`CG`).
* CIT v2 (`robust = TRUE`): test 3 is referred to a non-central F distribution whose
  non-centrality is estimated from the statistic for L in T ~ CT + G + L. A pair of
  variables can then not be called causal in both directions. Recommended when
  reverse causation is a concern; available for a single mediator. The default
  remains `robust = FALSE`.
* `test4 = c("empirical", "parametric")` selects how test 4's p-value is computed
  from null data. The default is now `"empirical"` for all functions (counting null
  statistics, up to `maxit = 10000`). `"parametric"` is the non-central F (or
  chi-square) method of Millstein et al. (2009), which `cit.cp()` used by default in
  2.3.2.

## Changes in results

* With a single instrument and a single mediator, `cit.cp()` and `cit.bp()` use the
  sign-pattern components by default (see above); set `sign.pattern = FALSE` for the
  2.3.2 components.

* `cit.cp()` now uses the empirical test 4 by default. Use `test4 = "parametric"` for
  the 2.3.2 method.
* The non-centrality in the parametric test 4 is the method-of-moments estimate
  lambda = F df1 (df2 - 2) / df2 - df1 (2.3.2 used df2 - df1 in place of df2 - 2;
  identical for two instruments, slightly different otherwise).
* The empirical test 4 p-value is (r + 1) / (B + 1) and cannot be 0.
* Each component test uses the observations complete for its own variables, as
  `cit.cp()` did in 2.3.2 (test 1: L, T, CT; test 2: L, G, T, CT; test 3: L, G, T,
  CG; test 4: all). With `robust = TRUE`, test 3 uses the test 4 observations.
* In permutation mode, each permutation is one reordering of all n observations,
  applied to every component test (L in tests 1 and 3 and in the v2 non-centrality,
  G in test 2, mediator-model residuals in test 4). In tests 1-3 the variable is
  permuted before incomplete observations are removed, so a shared `perm.index`
  gives the same permutation in every test and every CIT call.
* Permutation-mode test 3 in `cit.cp()` now permutes L (2.3.2 permuted G).
* Random numbers come from R, so `set.seed()` and `rseed` make every p-value
  reproducible; `rseed` no longer changes the caller's random number stream.
  (In 2.3.2, `cit.bp()` test 4 used a fixed internal seed and ignored `rseed`.)

## Bug fixes

* `cit.bp()` with permutations no longer crashes R when G or L has missing values,
  and no longer mixes the missing-value code into test 4's null data.
* Collinear instrument columns are dropped with correct degrees of freedom; rank
  decisions use a relative tolerance, so results do not depend on the scale of L.

## Internals

* GSL is no longer required. The test 4 resampling engine is a single C++ routine
  (RcppArmadillo); all other computation is in R.
* Tests (`tinytest`) compare component p-values with `lm()`, `glm()`, `coxph()` and
  `manova()`.

## Compatibility

* Argument order up to `rseed` is unchanged from 2.3.2. `C` is accepted with a
  deprecation warning. Output names are unchanged (`p_LindTgvnG`).
* **Breaking change:** all three functions return a data frame. Without
  permutations it has one row, as `cit.bp()` did in 2.3.2; `cit.cp()` previously
  returned a named vector. Code such as `r[["p_cit"]]` and `r$p_cit` works with
  every function; `r["p_cit"]` or `r[1]` on a `cit.cp()` result now gives a
  one-column data frame rather than a number.
* `fdr.cit()` accepts results from all three functions.

# cit 2.3.0 (2021-05-19)

* Fixed memory leaks with gsl (freeing memory).
* Fixed a bug in cit.bp() when using C. It was creating design matrices of
  the wrong size. That was causing memory leaks.
* Revamp for returning to CRAN.
* Adding a few tests for cit.bp().
