## Nested-model tests used by the CIT components. Each returns
## c(stat, df1, df2, p, sgn): an F test (linear), a likelihood-ratio chi-square
## (logistic, Cox; df2 = NA) or Wilks' lambda via Rao's F (multivariate).
## sgn is the sign of the coefficient of the added column when exactly one
## column is added (linear and logistic), otherwise NA.
## Degrees of freedom come from the ranks of the fitted designs, so collinear
## columns are dropped exactly as lm() and glm() drop them.

na_test <- function(df1 = NA_real_, df2 = NA_real_)
  c(stat = NA_real_, df1 = df1, df2 = df2, p = NA_real_, sgn = NA_real_)

## sign of the last coefficient when one column is added; NA otherwise
last_sign <- function(coef, df1) if (df1 == 1) unname(sign(coef[length(coef)])) else NA_real_

## linear outcome: F test for the columns of X1 not in X0
ftest_lin <- function(y, X0, X1) {
  q0 <- qr(X0, tol = 1e-7); q1 <- qr(X1, tol = 1e-7)
  df1 <- q1$rank - q0$rank; df2 <- length(y) - q1$rank
  if (df1 < 1 || df2 < 1) return(na_test(df1, df2))
  rss0 <- sum(qr.resid(q0, y)^2); rss1 <- sum(qr.resid(q1, y)^2)
  F <- (max(rss0 - rss1, 0) / df1) / (rss1 / df2)
  c(stat = F, df1 = df1, df2 = df2, p = pf(F, df1, df2, lower.tail = FALSE),
    sgn = last_sign(qr.coef(q1, y), df1))
}

## binary outcome: likelihood-ratio test from glm.fit deviances
lrt_logit <- function(y, X0, X1) {
  f0 <- suppressWarnings(glm.fit(X0, y, family = binomial()))
  f1 <- suppressWarnings(glm.fit(X1, y, family = binomial()))
  df1 <- f1$rank - f0$rank
  if (df1 < 1) return(na_test(df1))
  if (!f0$converged || !f1$converged)
    warning("logistic regression did not converge; the p-value may be unreliable", call. = FALSE)
  chi <- max(f0$deviance - f1$deviance, 0)
  c(stat = chi, df1 = df1, df2 = NA_real_, p = pchisq(chi, df1, lower.tail = FALSE),
    sgn = last_sign(f1$coefficients, df1))
}

## survival outcome: Cox partial-likelihood ratio test (Efron ties); X0 may have 0 columns
cox_loglik <- function(X, y) {
  fit <- survival::coxph.fit(X, y, strata = NULL, offset = NULL, init = NULL,
                             control = survival::coxph.control(), weights = NULL,
                             method = "efron", rownames = NULL)
  list(loglik = fit$loglik, rank = sum(!is.na(fit$coefficients)))
}
lrt_cox <- function(y, X0, X1) {
  f1 <- cox_loglik(X1, y)
  f0 <- if (ncol(X0) == 0) list(loglik = f1$loglik[1], rank = 0) else cox_loglik(X0, y)
  df1 <- f1$rank - f0$rank
  if (df1 < 1) return(na_test(df1))
  chi <- max(2 * (f1$loglik[length(f1$loglik)] - f0$loglik[length(f0$loglik)]), 0)
  c(stat = chi, df1 = df1, df2 = NA_real_, p = pchisq(chi, df1, lower.tail = FALSE), sgn = NA_real_)
}

## multivariate outcome Y (several mediators): Wilks' lambda with Rao's F,
## computed as in stats:::summary.manova
wilks_test <- function(Y, X0, X1) {
  q0 <- qr(X0, tol = 1e-7); q1 <- qr(X1, tol = 1e-7)
  q <- q1$rank - q0$rank; df.res <- nrow(Y) - q1$rank
  if (q < 1 || df.res < 1) return(na_test(q))
  E <- crossprod(qr.resid(q1, Y))
  H <- crossprod(qr.resid(q0, Y)) - E
  rss.qr <- qr(E)
  if (rss.qr$rank < ncol(Y)) return(na_test(q))
  eig <- Re(eigen(qr.coef(rss.qr, H), symmetric = FALSE, only.values = TRUE)$values)
  p <- length(eig)
  lam <- prod(1 / (1 + eig))
  tmp1 <- df.res - 0.5 * (p - q + 1)
  tmp2 <- (p * q - 2) / 4
  tmp3 <- p^2 + q^2 - 5
  tmp3 <- if (tmp3 > 0) sqrt(((p * q)^2 - 4) / tmp3) else 1
  F <- ((lam^(-1 / tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2)) / p / q
  df1 <- p * q; df2 <- tmp1 * tmp3 - 2 * tmp2
  c(stat = F, df1 = df1, df2 = df2, p = pf(F, df1, df2, lower.tail = FALSE), sgn = NA_real_)
}

## dispatch on outcome family for the outcome (T) models
otest <- function(family, y, X0, X1) {
  switch(family,
         linear   = ftest_lin(y, X0, X1),
         logistic = lrt_logit(y, X0, X1),
         cox      = lrt_cox(y, X0, X1))
}

## ---- non-central distributions --------------------------------------------
## Non-centrality from an observed statistic by the method of moments:
##   F:   E[F] = df2 (df1 + lambda) / (df1 (df2 - 2))  =>  lambda = F df1 (df2 - 2) / df2 - df1
##   chi: E[X] = df1 + lambda                           =>  lambda = X - df1
## truncated at 0.
ncp_mom <- function(stat, df1, df2 = NA) {
  lam <- if (is.na(df2)) stat - df1 else stat * df1 * (df2 - 2) / df2 - df1
  pmax(lam, 0)
}

## Patnaik's approximation, used when the non-centrality is too large for the
## exact algorithms: chi'(df, lam) ~ c * chi(h)
patnaik <- function(df, ncp) list(c = (df + 2 * ncp) / (df + ncp), h = (df + ncp)^2 / (df + 2 * ncp))
ncp_exact_max <- 1e4

## cdf of the non-central F (df2 finite) or chi-square (df2 = NA)
pnc <- function(q, df1, df2 = NA, ncp = 0, lower.tail = TRUE, log.p = FALSE) {
  if (is.na(df2)) {
    if (ncp == 0) return(pchisq(q, df1, lower.tail = lower.tail, log.p = log.p))
    if (ncp <= ncp_exact_max) return(pchisq(q, df1, ncp = ncp, lower.tail = lower.tail, log.p = log.p))
    a <- patnaik(df1, ncp)
    return(pchisq(q / a$c, a$h, lower.tail = lower.tail, log.p = log.p))
  }
  if (ncp == 0) return(pf(q, df1, df2, lower.tail = lower.tail, log.p = log.p))
  if (ncp <= ncp_exact_max) return(pf(q, df1, df2, ncp = ncp, lower.tail = lower.tail, log.p = log.p))
  a <- patnaik(df1, ncp)
  pf(q * df1 / (a$c * a$h), a$h, df2, lower.tail = lower.tail, log.p = log.p)
}

## probit of the cdf, computed in the tail that keeps precision
z_nc <- function(q, df1, df2, ncp) {
  lo <- pnc(q, df1, df2, ncp, lower.tail = TRUE, log.p = TRUE)
  up <- pnc(q, df1, df2, ncp, lower.tail = FALSE, log.p = TRUE)
  ifelse(lo < up, qnorm(lo, log.p = TRUE), -qnorm(up, log.p = TRUE))
}
