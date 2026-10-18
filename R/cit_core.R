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
## v2 (robust = TRUE, single mediator, linear or binary outcome): test 3 is
## referred to a non-central F whose non-centrality is the method-of-moments
## estimate from the statistic for L in T ~ CT + G + L (the observed test 4
## statistic), so that a pair cannot be called causal in both directions.

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
  if (family == "logistic" && !all(d$T %in% c(0, 1)))
    stop("'T' must be coded 0/1 for cit.bp.", call. = FALSE)
  if (family == "cox" && !all(d$E %in% c(0, 1)))
    stop("'E' must be coded 0/1 (1 = event).", call. = FALSE)
  if (d$n < ncol(d$L) + ncol(cbind(d$CT, d$CG, d$G)) + 4)
    stop("too few complete observations.", call. = FALSE)

  with_seed(rseed, {
    if (!is.null(perm.index)) {
      perm.index <- check_perm_index(perm.index, d$n.all)
      n.perm <- ncol(perm.index)
    } else if (n.perm > 0) {
      perm.index <- random_perms(d$n.all, n.perm)
    }
    perm <- if (n.perm > 0) restrict_perm(perm.index, d$keep) else matrix(integer(0), d$n, 0)
    run_cit(family, d, robust, test4, maxit, n.resampl, perm)
  })
}

run_cit <- function(family, d, robust, test4, maxit, n.resampl, perm) {
  n <- d$n; P <- ncol(perm)
  y <- if (family == "cox") survival::Surv(d$T, d$E) else d$T
  one <- if (family == "cox") NULL else rep(1, n)
  AT <- cbind(one, d$CT); if (is.null(AT)) AT <- matrix(0, n, 0)    # outcome-model base
  AG <- cbind(1, d$T, d$CG)                                           # test-3 base
  has3 <- family != "cox"

  comp <- function(Lx, Gx) {
    t1 <- otest(family, y, AT, cbind(AT, Lx))
    t2 <- otest(family, y, cbind(AT, d$L), cbind(AT, d$L, Gx))
    t3 <- if (!has3) NULL
          else if (ncol(d$G) > 1) wilks_test(d$G, AG, cbind(AG, Lx))
          else ftest_lin(drop(d$G), AG, cbind(AG, Lx))
    list(t1 = t1, t2 = t2, t3 = t3)
  }

  ## observed components and test-4 statistic
  s4 <- t4_setup(family, d)
  with_L <- function(Lx) {               # test-4 setup with instrument Lx
    sj <- s4
    sj$L <- Lx[, independent_cols(cbind(s4$A, d$G), Lx), drop = FALSE]
    sj
  }
  obs <- comp(d$L, d$G)
  stat4 <- t4_observed(s4)
  p3 <- if (has3) obs$t3[["p"]] else NULL
  if (robust) p3 <- p3_v2(obs$t3, stat4, t4_df(s4))

  ## permuted components: L permuted for tests 1 and 3, G for test 2;
  ## for v2 the non-centrality comes from the same permuted L
  pc <- lapply(seq_len(P), function(j) {
    pj <- perm[, j]
    cj <- comp(d$L[pj, , drop = FALSE], d$G[pj, , drop = FALSE])
    p3j <- if (has3) cj$t3[["p"]] else NULL
    if (robust) {
      sj <- with_L(d$L[pj, , drop = FALSE])
      p3j <- p3_v2(cj$t3, t4_observed(sj), t4_df(sj))
    }
    c(cj$t1[["p"]], cj$t2[["p"]], p3j)
  })

  ## test 4
  maxp <- max(obs$t1[["p"]], obs$t2[["p"]], p3)
  if (test4 == "empirical") {
    obs_perm <- if (P > 0) t4_stats(s4, perm) else numeric(0)
    r4 <- t4_empirical(s4, stat4, obs_perm, maxp, maxit)
  } else {
    r4 <- t4_parametric(s4, stat4, perm, n.resampl)
  }

  out <- c(obs$t1[["p"]], obs$t2[["p"]], p3, r4$p)
  nm <- if (has3) pnames[-1] else pnames[c(2, 3, 5)]
  res <- c(max(out), out); names(res) <- c("p_cit", nm)
  if (P == 0) return(res)

  pm <- cbind(do.call(rbind, pc), r4$p_perm)
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
