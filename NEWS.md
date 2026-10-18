# cit 3.0.0

## New features

* `cit.cox()` for a censored survival outcome. The test of G given T is omitted:
  the mediator must be measured before the events occur, which rules out reverse
  causation by design (see `?cit.cox`).
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

* `cit.cp()` now uses the empirical test 4 by default. Use `test4 = "parametric"` for
  the 2.3.2 method.
* The non-centrality in the parametric test 4 is the method-of-moments estimate
  lambda = F df1 (df2 - 2) / df2 - df1 (2.3.2 used df2 - df1 in place of df2 - 2;
  identical for two instruments, slightly different otherwise).
* The empirical test 4 p-value is (r + 1) / (B + 1) and cannot be 0.
* Observations with any missing value are removed once, so all component tests use
  the same observations (2.3.2 `cit.cp()` used the available cases for each test).
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
  deprecation warning. Output names are unchanged (`p_LindTgvnG`); without
  permutations `cit.cp()` and `cit.cox()` return a named vector and `cit.bp()` a
  one-row data frame, as in 2.3.2.
* `fdr.cit()` accepts results from all three functions.
