## Sign check as its own column, sign_ind (3.1.0)

set.seed(1); n <- 423
L <- rbinom(n, 2, 0.3); U <- rnorm(n)

## signs agree, (-, +, -), with a strong effect of L on T: the robust test 3 p-value is
## numerically 1. In 3.0.0 the sign-pattern tie-break then also set test 1 to 1.
A <- -1.2 * L + 0.5 * U + rnorm(n)
B <- -0.4 * L + 0.5 * U + rnorm(n)
r  <- cit.cp(L, B, A, robust = TRUE, rseed = 1)
r0 <- cit.cp(L, B, A, robust = TRUE, sign.pattern = FALSE, rseed = 1)
expect_equal(r$sign_ind, 0)
expect_equal(r0$p_GassocLgvnT, 1)
expect_equal(unlist(r[2:5]), unlist(r0[2:5]))          # robust: two-sided components
expect_true(r$p_TassocL < 1e-30)
expect_equal(r$p_cit, max(unlist(r0[2:5])))
expect_false("sign_ind" %in% names(r0))

## signs agree, robust = FALSE: tests 1-3 report half their two-sided p-values
s  <- cit.cp(L, B, A, rseed = 1, test4 = "parametric")
s0 <- cit.cp(L, B, A, sign.pattern = FALSE, rseed = 1, test4 = "parametric")
expect_equal(s$sign_ind, 0)
expect_equal(unlist(s[2:4]), unlist(s0[2:4]) / 2)
expect_equal(s$p_LindTgvnG, s0$p_LindTgvnG)
expect_equal(s$p_cit, max(unlist(s[2:5])))

## signs disagree, (+, -, +): a confounder U of G and T with opposite loadings
G <- 0.8 * L + U + rnorm(n)
T <- 0.5 * L - U + rnorm(n)
v  <- cit.cp(L, G, T, rseed = 1, test4 = "parametric")
v0 <- cit.cp(L, G, T, sign.pattern = FALSE, rseed = 1, test4 = "parametric")
expect_equal(v$sign_ind, 1)
expect_equal(v$p_cit, 1)
expect_equal(unlist(v[2:4]), unlist(v0[2:4]) / 2)      # components are not altered
expect_true(all(unlist(v[2:4]) < 1))
vr <- cit.cp(L, G, T, robust = TRUE, rseed = 1, test4 = "parametric")
expect_equal(vr$sign_ind, 1)
expect_equal(vr$p_cit, 1)
expect_true(all(unlist(vr[2:3]) < 1))

## binary outcome
vb <- cit.bp(L, G, as.numeric(T > median(T)), rseed = 1, test4 = "parametric")
expect_equal(vb$sign_ind, 1)
expect_equal(vb$p_cit, 1)

## several mediators: signs not defined, no sign_ind
m <- cit.cp(L, cbind(G, B), T, rseed = 1, test4 = "parametric")
expect_false("sign_ind" %in% names(m))
expect_equal(m$p_cit, max(unlist(m[2:5])))

## sign_check: an undefined sign gives the two-sided p-values and NA
expect_equal(cit:::sign_check(c(0.2, 0.1, 0.4), c(1, 0, 1)), c(0.2, 0.1, 0.4, NA))
expect_equal(cit:::sign_check(c(0.2, 0.1, 0.4), c(1, -1, -1)), c(0.1, 0.05, 0.2, 0))
expect_equal(cit:::sign_check(c(0.2, 0.1, 0.4), c(-1, -1, -1), halve = FALSE), c(0.2, 0.1, 0.4, 1))

## permutation rows carry sign_ind, and p_cit includes it
rp <- cit.cp(L, G, T, n.perm = 10, rseed = 1, test4 = "parametric")
expect_equal(names(rp), c("perm", "p_cit", "p_TassocL", "p_TassocGgvnL", "p_GassocLgvnT",
                          "p_LindTgvnG", "sign_ind"))
expect_true(all(rp$sign_ind %in% c(0, 1)))
expect_equal(rp$p_cit, apply(rp[, 3:7], 1, max))

## cit.cox has no sign check
cx <- cit.cox(L, G, rexp(n), rbinom(n, 1, 0.7), rseed = 1, test4 = "parametric")
expect_false("sign_ind" %in% names(cx))

## fdr.cit: sign_ind applies to the observed tests only
pidx <- replicate(20, sample.int(n))
res <- lapply(1:12, function(i) {
  set.seed(100 + i); U <- rnorm(n); G <- 0.6 * L + U + rnorm(n)
  T <- if (i %% 3 == 0) 0.5 * L - U + rnorm(n) else 0.5 * G + rnorm(n)
  cit.cp(L, G, T, perm.index = pidx, rseed = i, test4 = "parametric")
})
bad <- vapply(res, function(x) x$sign_ind[1] == 1, NA)
f  <- fdr.cit(res)
fn <- fdr.cit(lapply(res, function(x) x[, names(x) != "sign_ind"]))
expect_true(any(bad) && any(!bad))
expect_equal(f$sign_ind, as.numeric(bad))
expect_true(all(f$q.cit[bad] == 1 & f$q.cit.ll[bad] == 1 & f$q.cit.ul[bad] == 1))
expect_true(all(f$p.cit.cal[bad] == 1))
expect_equal(f$q.cit[!bad], fn$q.cit[!bad])            # the others: as without the sign check
expect_equal(f$p.cit.cal[!bad], fn$p.cit.cal[!bad])
expect_equal(f$q.TaL, fn$q.TaL)                        # component q-values ignore sign_ind
expect_false("sign_ind" %in% names(fn))
