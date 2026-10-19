## One implementation of the causal inference test for every outcome type.
##
## Components (L instrument, G mediator, T outcome; CT covariates for T models,
## CG covariates for G models):
##   1. T ~ L | CT                 association of T with L
##   2. T ~ G | L, CT              association of T with G given L
##   3. G ~ L | T, CG              association of G with L given T   (not for Cox)
##   4. L ind T | G                see test4.R
## p_cit is the maximum of the component p-values (intersection-union test).
##
## Missing data: each test uses the rows complete for its own variables
## (test 1: L, T, CT; test 2: L, G, T, CT; test 3: L, G, T, CG; test 4: all).
##
## v2 (robust = TRUE, single mediator, linear or binary outcome): test 3 is
## referred to a non-central F whose non-centrality is the method-of-moments
## estimate from the statistic for L in T ~ CT + G + L (the observed test 4
## statistic), so that a pair cannot be called causal in both directions.
## Test 3 then uses the test 4 rows, so both statistics come from the same data.
##
## Permutations: permutation j is one reordering of all n rows, applied to L in
## tests 1 and 3 (and in the v2 non-centrality), to G in test 2 and to the
## mediator-model residuals in test 4, so dependence among tests, and among CITs
## sharing a perm.index, is preserved. In tests 1-3 the variable is permuted
## first (missing values travel with their rows) and incomplete rows are dropped
## afterwards; test 4 restricts the permutation to its rows (restrict_perm).

pnames <- c("p_cit", "p_TassocL", "p_TassocGgvnL", "p_GassocLgvnT", "p_LindTgvnG")

cit_core <- function(family, L, G, T, E = NULL, CT = NULL, CG = NULL, robust = FALSE,
                     test4 = "empirical", maxit = 10000, n.resampl = 50, n.perm = 0,
                     perm.index = NULL, rseed = NULL) {
  d <- prep_data(L, G, T, CT, CG, E)
  if (is.null(d$L) || ncol(d$L) < 1) stop("'L' must have at least one column.", call. = FALSE)
  m <- ncol(d$G)
  if (robust && m > 1)
    stop("robust = TRUE (CIT v2) is not yet available for multiple mediators; use robust = FALSE.", call. = FALSE)
  if (robust && family == "cox")
    stop("robust = TRUE does not apply to cit.cox, which has no test 3.", call. = FALSE)
  if (family == "logistic" && !all(d$T %in% c(0, 1, NA)))
    stop("'T' must be coded 0/1 for cit.bp.", call. = FALSE)
  if (family == "cox" && !all(d$E %in% c(0, 1, NA)))
    stop("'E' must be coded 0/1 (1 = event).", call. = FALSE)
  if (sum(Reduce(`&`, d$ok)) < ncol(d$L) + ncol(cbind(d$CT, d$CG, d$G)) + 4)
    stop("too few complete observations.", call. = FALSE)

  with_seed(rseed, {
    if (!is.null(perm.index)) {
      perm <- check_perm_index(perm.index, d$n.all)
    } else {
      perm <- random_perms(d$n.all, n.perm)
    }
    run_cit(family, d, robust, test4, maxit, n.resampl, perm)
  })
}

run_cit <- function(family, d, robust, test4, maxit, n.resampl, perm) {
  P <- ncol(perm); ok <- d$ok
  has3 <- family != "cox"
  sub <- function(x, r) if (is.null(x)) NULL else x[r, , drop = FALSE]
  yv <- function(r) if (family == "cox") survival::Surv(d$T[r], d$E[r]) else d$T[r]
  AT <- function(r) {                                   # outcome-model base
    A <- cbind(if (family != "cox") rep(1, sum(r)), sub(d$CT, r))
    if (is.null(A)) matrix(0, sum(r), 0) else A
  }
  r4 <- Reduce(`&`, ok)                                 # rows for test 4

  ## components 1-3 with instrument Lx and mediator Gx in test 2 (permuted or not)
  comp <- function(Lx, okL, Gx, okG) {
    r1 <- okL & ok$T & ok$CT
    A1 <- AT(r1)
    t1 <- otest(family, yv(r1), A1, cbind(A1, sub(Lx, r1)))
    r2 <- ok$L & okG & ok$T & ok$CT
    A2 <- cbind(AT(r2), sub(d$L, r2))
    t2 <- otest(family, yv(r2), A2, cbind(A2, sub(Gx, r2)))
    if (!has3) return(c(t1[["p"]], t2[["p"]]))
    r3 <- okL & ok$G & ok$T & ok$CG
    if (robust) r3 <- r3 & ok$CT
    AG <- cbind(1, d$T[r3], sub(d$CG, r3))
    G3 <- sub(d$G, r3); X3 <- cbind(AG, sub(Lx, r3))
    t3 <- if (ncol(G3) > 1) wilks_test(G3, AG, X3) else ftest_lin(drop(G3), AG, X3)
    p3 <- t3[["p"]]
    if (robust) {                     # non-centrality from the same rows and the same Lx
      s <- t4_setup(family, sub_data(d, r3, Lx))
      p3 <- p3_v2(t3, t4_observed(s), t4_df(s))
    }
    c(t1[["p"]], t2[["p"]], p3)
  }

  obs <- comp(d$L, ok$L, d$G, ok$G)
  pc <- lapply(seq_len(P), function(j) {
    pj <- perm[, j]
    comp(d$L[pj, , drop = FALSE], ok$L[pj], d$G[pj, , drop = FALSE], ok$G[pj])
  })

  ## test 4
  s4 <- t4_setup(family, sub_data(d, r4))
  stat4 <- t4_observed(s4)
  perm4 <- if (P > 0) restrict_perm(perm, r4) else matrix(integer(0), s4$n, 0)
  if (test4 == "empirical") {
    obs_perm <- if (P > 0) t4_stats(s4, perm4) else numeric(0)
    res4 <- t4_empirical(s4, stat4, obs_perm, max(obs), maxit)
  } else {
    res4 <- t4_parametric(s4, stat4, perm4, n.resampl)
  }

  out <- c(obs, res4$p)
  nm <- if (has3) pnames[-1] else pnames[c(2, 3, 5)]
  res <- c(max(out), out); names(res) <- c("p_cit", nm)
  if (P == 0) return(res)

  pm <- cbind(do.call(rbind, pc), res4$p_perm)
  pm <- cbind(apply(pm, 1, max), pm)
  df <- as.data.frame(rbind(res, pm), row.names = FALSE)
  names(df) <- names(res)
  cbind(perm = 0:P, df)
}

## v2 test 3: F for G ~ L | T, CG against a non-central F with lambda from the
## test 4 statistic (F, or likelihood-ratio chi-square for a binary outcome)
p3_v2 <- function(t3, stat4, df4) {
  if (anyNA(t3[c("stat", "df1", "df2")]) || is.na(stat4)) return(NA_real_)
  lam <- ncp_mom(stat4, df4[["df1"]], df4[["df2"]])
  pnc(t3[["stat"]], t3[["df1"]], t3[["df2"]], lam, lower.tail = FALSE)
}
