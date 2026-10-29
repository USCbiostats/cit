## Component test 4: is the association of L with T given G weaker than expected
## if G were not a mediator?
##
## Null data are generated as G* = fitted(G ~ L + CG) + row-permuted residuals:
## G* keeps its association with L but has no link to T other than through L,
## so in T ~ CT + G* + L the whole effect of L remains. The observed statistic
## for L in T ~ CT + G + L is compared with the statistics from G*. Small values
## (L explained by G) give small p-values.
##
## method = "empirical": p4 = (#{stat* <= stat} + 1) / (B + 1), with B up to
##   maxit; resampling stops early unless p4 could decide p_cit (rule from cit.bp
##   in cit <= 2.3.2: 1000 first, continue while fewer than 20 null statistics fall
##   at or below the observed one, the other components are < 0.1 and smaller than
##   the current p4 bound).
## method = "parametric" (Millstein et al. 2009): n.resampl null statistics give
##   the non-centrality by the method of moments; statistics are mapped to z
##   through the fitted non-central cdf, and p4 = Phi((z - mean(z*)) / sd(z*)).

firstloop <- 1000L
posno <- 20L
alpha_ext <- 0.1

## Everything test 4 needs, fixed before resampling
t4_setup <- function(family, d) {
  n <- d$n
  A <- if (family == "cox") d$CT else cbind(rep(1, n), d$CT)
  if (is.null(A)) A <- matrix(0, n, 0)
  # keep linearly independent columns only (as lm() does for tests 1-3): A, then
  # G given A, then L given A and G
  A <- A[, independent_cols(matrix(0, n, 0), A), drop = FALSE]
  G <- d$G[, independent_cols(A, d$G), drop = FALSE]
  BG <- cbind(1, d$L, d$CG)                         # G model: G ~ L + CG
  qg <- qr(BG, tol = 1e-7)
  Ghat <- qr.fitted(qg, G); Gres <- qr.resid(qg, G)
  Lk <- d$L[, independent_cols(cbind(A, G), d$L), drop = FALSE]
  list(family = family, y = d$T, E = d$E, A = A, L = Lk,
       Ghat = as.matrix(Ghat), Gres = as.matrix(Gres), n = n)
}

## statistics for the columns of `perm` (n x B, 1-based row indices)
t4_stats <- function(s, perm) {
  if (ncol(s$L) == 0) return(rep(NA_real_, ncol(perm)))   # L not estimable given CT and G
  if (s$family == "cox") {
    y <- survival::Surv(s$y, s$E)
    return(apply(perm, 2, function(p) {
      Gs <- s$Ghat + s$Gres[p, , drop = FALSE]
      X0 <- cbind(s$A, Gs)
      tryCatch(suppressWarnings(lrt_cox(y, X0, cbind(X0, s$L))[["stat"]]),
               error = function(e) NA_real_)
    }))
  }
  .boot_test4(s$y, s$A, s$L, s$Ghat, s$Gres, perm - 1L,
              family = if (s$family == "linear") 0L else 1L)
}

t4_df <- function(s) {
  df1 <- ncol(s$L)
  df2 <- if (s$family == "linear") s$n - ncol(s$A) - ncol(s$Ghat) - df1 else NA_real_
  c(df1 = df1, df2 = df2)
}

## observed statistic: G* = G (identity permutation), same code path as the null
t4_observed <- function(s) t4_stats(s, matrix(seq_len(s$n), ncol = 1))

random_perms <- function(n, B) {
  if (B < 1) return(matrix(integer(0), n, 0))
  vapply(seq_len(B), function(i) sample.int(n), integer(n))
}

## empirical p-values for the observed statistic and the permuted statistics
## `obs_perm` (one per perm.index column); `maxp` = max of the other components
t4_empirical <- function(s, obs, obs_perm, maxp, maxit) {
  B <- min(firstloop, maxit)
  st <- t4_stats(s, random_perms(s$n, B))
  le <- function(stat, o) sum(is.na(stat) | stat <= o)   # failed fits count toward the null
  npos <- le(st, obs)
  extend <- !is.na(obs) && npos < posno && !is.na(maxp) && maxp < alpha_ext &&
            maxp < (npos + 1) / B && B < maxit
  while (extend) {
    more <- t4_stats(s, random_perms(s$n, min(firstloop, maxit - length(st))))
    hits <- cumsum(is.na(more) | more <= obs)
    stop_at <- which(npos + hits >= posno)[1]
    if (!is.na(stop_at)) more <- more[seq_len(stop_at)]
    st <- c(st, more)
    npos <- le(st, obs)
    extend <- npos < posno && length(st) < maxit
  }
  B <- length(st)
  p <- if (is.na(obs)) NA_real_ else (npos + 1) / (B + 1)
  pp <- vapply(obs_perm, function(o) if (is.na(o)) NA_real_ else (le(st, o) + 1) / (B + 1), 0)
  list(p = p, p_perm = pp, B = B)
}

## parametric p-value of `obs` given null statistics `st`. A statistic at the edge of the
## support (e.g. exactly 0, after rounding, when the fitted non-centrality is 0) maps to an
## infinite z and is left out of the mean and sd; p is kept above 0.
t4_para_p <- function(obs, st, df) {
  st <- st[is.finite(st)]
  if (is.na(obs) || length(st) < 3) return(NA_real_)
  lam <- ncp_mom(mean(st), df[["df1"]], df[["df2"]])
  z <- z_nc(st, df[["df1"]], df[["df2"]], lam)
  z <- z[is.finite(z)]
  if (length(z) < 3) return(NA_real_)
  z0 <- z_nc(obs, df[["df1"]], df[["df2"]], lam)
  max(pnorm(z0, mean = mean(z), sd = sd(z)), .Machine$double.xmin)
}

## parametric p-values; the first n.perm null statistics use perm.index, and the
## permuted p-value for column j treats statistic j as observed (leave-one-out)
t4_parametric <- function(s, obs, perm, n.resampl) {
  P <- ncol(perm)
  R <- max(n.resampl, P, if (P > 0) 4L else 3L)   # leave-one-out keeps >= 3
  st <- t4_stats(s, cbind(perm, random_perms(s$n, R - P)))
  df <- t4_df(s)
  p <- t4_para_p(obs, st, df)
  pp <- vapply(seq_len(P), function(j) t4_para_p(st[j], st[-j], df), 0)
  list(p = p, p_perm = pp, B = R, obs_perm = st[seq_len(P)])
}
