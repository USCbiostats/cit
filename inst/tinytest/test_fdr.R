## fdr.cit: omnibus q-value at the common threshold p_cit (3.0.0)
set.seed(20261021)
m <- 300; K <- 4; B <- 10
pn <- c("p_TassocL", "p_TassocGgvnL", "p_GassocLgvnT", "p_LindTgvnG")
H <- matrix(rbinom(m * K, 1, 0.5), m, K)                       # 1 = component null
P <- ifelse(H == 1, runif(m * K), pnorm(rnorm(m * K) + 3, lower.tail = FALSE))
Pp <- lapply(seq_len(B), function(b) { d <- as.data.frame(matrix(runif(m * K), m, K)); names(d) <- pn; d })
mk <- function(P, Pp, pn) lapply(seq_len(nrow(P)), function(i) {
  x <- rbind(P[i, ], t(vapply(Pp, function(d) unlist(d[i, pn]), numeric(length(pn)))))
  d <- as.data.frame(x); names(d) <- pn
  cbind(p_cit = apply(x, 1, max), d)
})
f <- fdr.cit(mk(P, Pp, pn))

## shape and bounds
expect_equal(nrow(f), m)
expect_true(all(c("p.cit", "q.cit", "q.cit.ll", "q.cit.ul", "q.TaL", "q.LiTgvG") %in% names(f)))
expect_true(all(f$q.cit >= 0 & f$q.cit <= 1))
expect_true(all(f$q.cit.ll <= f$q.cit + 1e-12 & f$q.cit <= f$q.cit.ul + 1e-12))

## q.cit is non-decreasing in p.cit
o <- order(f$p.cit)
expect_true(all(diff(f$q.cit[o]) >= -1e-12))

## q.cit equals a direct implementation: component FDRs at t = p_cit, product, running minimum
Pd <- as.data.frame(P); names(Pd) <- pn
pc <- apply(P, 1, max)
Q <- vapply(pc, function(t) 1 - prod(vapply(pn, function(k) 1 - min(1, cit:::fdr.od(Pd[[k]], Pp, k, m, t)[1]), 0)), 0)
oo <- order(pc, decreasing = TRUE); Q[oo] <- cummin(Q[oo])
expect_equal(f$q.cit, Q, tolerance = 1e-10)

## never smaller than the 2.3.2 combination of the (monotone) component q-values, before its running minimum
q_own <- 1 - apply(1 - as.matrix(f[, c("q.TaL", "q.TaGgvL", "q.GaLgvT", "q.LiTgvG")]), 1, prod)
q_own[order(pc)] <- rev(cummin(rev(q_own[order(pc)])))
expect_true(mean(f$q.cit >= q_own - 1e-10) > 0.99)

## three components (cit.cox with test3 = FALSE)
pn3 <- pn[-3]
f3 <- fdr.cit(mk(P[, -3], lapply(Pp, `[`, pn3), pn3))
expect_false("q.GaLgvT" %in% names(f3))
expect_true(all(diff(f3$q.cit[order(f3$p.cit)]) >= -1e-12))
