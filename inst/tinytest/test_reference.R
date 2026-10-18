# Component p-values against lm(), glm(), coxph() and manova() on simulated data

sim <- function(seed, n = 200, family = "linear", m = 1) {
  set.seed(seed)
  L <- matrix(rbinom(n * 2, 2, .5), n)
  CT <- matrix(rnorm(n * 2), n); CG <- matrix(rnorm(n), n)
  G1 <- drop(L %*% c(.4, .3)) + .3 * CG[, 1] + rnorm(n)
  G <- if (m == 1) matrix(G1) else cbind(G1, .5 * G1 + rnorm(n))
  eta <- .4 * G1 + drop(CT %*% c(.3, .3))
  out <- list(L = L, G = G, CT = CT, CG = CG)
  if (family == "linear") out$T <- eta + rnorm(n)
  if (family == "logistic") out$T <- rbinom(n, 1, plogis(eta))
  if (family == "cox") { out$T <- rexp(n, exp(eta)); out$E <- rbinom(n, 1, .8) }
  out
}

# continuous outcome, with a missing value in T
d <- sim(1); d$T[3] <- NA
r <- cit.cp(d$L, d$G, d$T, d$CT, CG = d$CG, rseed = 1, maxit = 200)
k <- !is.na(d$T)
L <- d$L[k, ]; G <- d$G[k, 1]; T <- d$T[k]; CT <- d$CT[k, ]; CG <- d$CG[k, ]
expect_equal(r[["p_TassocL"]], anova(lm(T ~ CT), lm(T ~ CT + L))$`Pr(>F)`[2])
expect_equal(r[["p_TassocGgvnL"]], anova(lm(T ~ CT + L), lm(T ~ CT + L + G))$`Pr(>F)`[2])
expect_equal(r[["p_GassocLgvnT"]], anova(lm(G ~ T + CG), lm(G ~ T + CG + L))$`Pr(>F)`[2])

# v2: non-central F with method-of-moments lambda
r2 <- cit.cp(d$L, d$G, d$T, d$CT, CG = d$CG, rseed = 1, robust = TRUE, maxit = 200)
a5 <- anova(lm(T ~ CT + G), lm(T ~ CT + G + L))
lam <- max(0, a5$F[2] * a5$Df[2] * (a5$Res.Df[2] - 2) / a5$Res.Df[2] - a5$Df[2])
a3 <- anova(lm(G ~ T + CG), lm(G ~ T + CG + L))
expect_equal(r2[["p_GassocLgvnT"]], pf(a3$F[2], a3$Df[2], a3$Res.Df[2], ncp = lam, lower.tail = FALSE))

# binary outcome: tests 1-2 are likelihood-ratio tests, test 3 is linear
d <- sim(2, family = "logistic")
r <- cit.bp(d$L, d$G, d$T, d$CT, rseed = 1, maxit = 200)
L <- d$L; G <- d$G[, 1]; T <- d$T; CT <- d$CT
expect_equal(r[["p_TassocL"]], anova(glm(T ~ CT, family = binomial), glm(T ~ CT + L, family = binomial), test = "LRT")$`Pr(>Chi)`[2])
expect_equal(r[["p_GassocLgvnT"]], anova(lm(G ~ T), lm(G ~ T + L))$`Pr(>F)`[2])

# survival outcome: three components
d <- sim(3, family = "cox")
r <- cit.cox(d$L, d$G, d$T, d$E, d$CT, rseed = 1, maxit = 200)
expect_equal(names(r), c("p_cit", "p_TassocL", "p_TassocGgvnL", "p_LindTgvnG"))
f0 <- survival::coxph(survival::Surv(d$T, d$E) ~ d$CT + d$L)
f1 <- survival::coxph(survival::Surv(d$T, d$E) ~ d$CT + d$L + d$G)
expect_equal(r[["p_TassocGgvnL"]], pchisq(2 * (f1$loglik[2] - f0$loglik[2]), 1, lower.tail = FALSE), tolerance = 1e-6)

# multiple mediators: test 3 is Wilks' lambda
d <- sim(4, m = 2)
r <- cit.cp(d$L, d$G, d$T, rseed = 1, maxit = 200)
Lm <- d$L; X <- d$T
expect_equal(r[["p_GassocLgvnT"]], summary(manova(d$G ~ X + Lm), test = "Wilks")$stats["Lm", "Pr(>F)"])

# missing values in permutation mode do not crash and give complete results
d <- sim(5, m = 2); d$G[7, 2] <- NA; d$L[9, 1] <- NA
r <- cit.cp(d$L, d$G, d$T, d$CT, n.perm = 3, rseed = 1, maxit = 200)
expect_equal(dim(r), c(4L, 6L))
expect_false(anyNA(r))

# reproducibility
d <- sim(6)
expect_identical(cit.cp(d$L, d$G, d$T, rseed = 3), cit.cp(d$L, d$G, d$T, rseed = 3))
set.seed(4); a <- cit.cp(d$L, d$G, d$T); set.seed(4); b <- cit.cp(d$L, d$G, d$T)
expect_identical(a, b)

# p4 is never exactly 0
expect_true(cit.cp(d$L, d$G, d$T, rseed = 1, maxit = 100)[["p_LindTgvnG"]] > 0)
