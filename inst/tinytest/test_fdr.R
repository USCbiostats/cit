## fdr.cit: omnibus q-value on the calibrated scale (3.0.0)
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
expect_true(all(c("p.cit", "p.cit.cal", "q.cit", "q.cit.ll", "q.cit.ul", "q.TaL", "q.LiTgvG") %in% names(f)))
expect_true(all(f$q.cit >= 0 & f$q.cit <= 1))
expect_true(all(f$q.cit.ll <= f$q.cit + 1e-12 & f$q.cit <= f$q.cit.ul + 1e-12))

## q.cit is non-decreasing in p.cit.cal
o <- order(f$p.cit.cal)
expect_true(all(diff(f$q.cit[o]) >= -1e-12))

## p.cit.cal and q.cit equal a direct implementation: each component mapped to its pooled permutation
## null CDF, component FDRs over the tests with F_k(p_k) <= u, product, running minimum
Pd <- as.data.frame(P); names(Pd) <- pn
Fk <- sapply(pn, function(k) ecdf(unlist(lapply(Pp, `[[`, k)))(Pd[[k]]))
u <- apply(Fk, 1, max)
expect_equal(f$p.cit.cal, u, tolerance = 1e-12)
Q <- vapply(u, function(v) 1 - prod(vapply(seq_along(pn), function(k)
  1 - min(1, cit:::fdr.od(Pd[[k]], Pp, pn[k], m, max(Pd[[k]][Fk[, k] <= v]))[1]), 0)), 0)
oo <- order(u, decreasing = TRUE); Q[oo] <- cummin(Q[oo])
expect_equal(f$q.cit, Q, tolerance = 1e-10)

## component nulls on different ranges (test 1 null U(0, .49), test 2 null U(.5, 1)): calibration
## puts both on the same null scale, so tests non-null for both components are called
m2 <- 600; pn2 <- pn[1:2]
H2 <- matrix(rbinom(m2 * 2, 1, 0.5), m2, 2)
sc <- function(p, k) if (k == 1) 0.49 * p else 0.5 + 0.5 * p
P2 <- sapply(1:2, function(k) sc(ifelse(H2[, k] == 1, runif(m2), pnorm(rnorm(m2) + 3, lower.tail = FALSE)), k))
Pp2 <- lapply(seq_len(B), function(b) { d <- as.data.frame(sapply(1:2, function(k) sc(runif(m2), k))); names(d) <- pn2; d })
f2 <- fdr.cit(mk(P2, Pp2, pn2))
alt2 <- H2[, 1] == 0 & H2[, 2] == 0
expect_true(median(f2$q.cit[alt2]) < 0.1)
expect_true(mean(!alt2[f2$q.cit <= 0.1]) < 0.2)

## three components (cit.cox with test3 = FALSE)
pn3 <- pn[-3]
f3 <- fdr.cit(mk(P[, -3], lapply(Pp, `[`, pn3), pn3))
expect_false("q.GaLgvT" %in% names(f3))
expect_true(all(diff(f3$q.cit[order(f3$p.cit.cal)]) >= -1e-12))
