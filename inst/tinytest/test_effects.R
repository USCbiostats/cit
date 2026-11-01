## cit.effects: estimates against lm/glm/coxph fits, and input handling
library(survival)
set.seed(21); n <- 800
L <- rbinom(n, 2, .3); C1 <- rnorm(n); G <- .4 * L + .3 * C1 + rnorm(n)
Tc <- .3 * G + .1 * L + .2 * C1 + rnorm(n)
Tb <- rbinom(n, 1, plogis(-3 + .4 * G + .2 * L))
Ty <- rexp(n, .02 * exp(.4 * G + .2 * L)); Ev <- rbinom(n, 1, .2)

## continuous: product of coefficients and Sobel-type SE
e <- cit.effects(L, G, Tc, CT = C1, CG = C1)
fm <- lm(G ~ L + C1); fo <- lm(Tc ~ L + G + C1)
b1 <- coef(fm)[["L"]]; t1 <- coef(fo)[["L"]]; t2 <- coef(fo)[["G"]]
expect_equal(e$estimate[e$effect == "direct"], t1)
expect_equal(e$estimate[e$effect == "indirect"], t2 * b1)
expect_equal(e$estimate[e$effect == "total"], t1 + t2 * b1)
expect_equal(e$estimate[e$effect == "prop.mediated"], t2 * b1 / (t1 + t2 * b1))
expect_equal(e$se[e$effect == "indirect"], sqrt(t2^2 * vcov(fm)["L", "L"] + b1^2 * vcov(fo)["G", "G"]))
expect_equal(e$scale[1:3], rep("difference", 3))
expect_true(all(e$lower < e$estimate & e$estimate < e$upper))

## binary: odds ratios
e <- cit.effects(L, G, Tb, outcome = "binary")
fo <- glm(Tb ~ L + G, family = binomial); b1 <- coef(lm(G ~ L))[["L"]]
expect_equal(e$estimate[e$effect == "direct"], exp(coef(fo)[["L"]]))
expect_equal(e$estimate[e$effect == "indirect"], exp(coef(fo)[["G"]] * b1))
expect_equal(e$se[e$effect == "direct"], sqrt(vcov(fo)["L", "L"]))
expect_equal(e$scale[1], "OR")

## survival: hazard ratios (Efron ties), rare-outcome warning
expect_warning(e <- cit.effects(L, G, Ty, E = Ev, outcome = "survival"), "rare outcome")
fo <- coxph(Surv(Ty, Ev) ~ L + G); b1 <- coef(lm(G ~ L))[["L"]]
expect_equal(e$estimate[e$effect == "direct"], exp(unname(coef(fo)["L"])))
expect_equal(e$estimate[e$effect == "indirect"], exp(unname(coef(fo)["G"]) * b1))
expect_equal(e$scale[1], "HR")

## several columns of L: one block of four rows per column, named
L2 <- cbind(snp1 = L, snp2 = rbinom(n, 2, .4))
e <- cit.effects(L2, G, Tc)
expect_equal(nrow(e), 8)
expect_equal(unique(e$L), c("snp1", "snp2"))
expect_equal(e$estimate[e$L == "snp2" & e$effect == "direct"], coef(lm(Tc ~ L2 + G))[["L2snp2"]])

## missing values: complete cases across all variables
Gm <- G; Gm[1:10] <- NA; Tm <- Tc; Tm[11:20] <- NA
e <- cit.effects(L, Gm, Tm)
expect_equal(e$n[1], n - 20)
ok <- !is.na(Gm) & !is.na(Tm)
expect_equal(e$estimate[e$effect == "direct"], coef(lm(Tm[ok] ~ L[ok] + Gm[ok]))[[2]])

## input checks
expect_error(cit.effects(L, cbind(G, G), Tc), "single mediator")
expect_error(cit.effects(L, G, Ty, outcome = "survival"), "event indicator")
expect_error(cit.effects(L, G, Tc, E = Ev), "only with")
expect_error(cit.effects(L, G, Tc, outcome = "binary"), "0/1")
expect_error(cit.effects(L, G, Tc, conf.level = 1.5), "conf.level")
