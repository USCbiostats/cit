## Direct and indirect effect estimates to accompany the CIT, by the regression
## (product-of-coefficients) estimators of VanderWeele and Vansteelandt, in the
## form given by Valeri and VanderWeele (2013) for models without an
## exposure-mediator interaction (the CIT's models have none). Standard errors
## are by the delta method, treating the mediator and outcome models as
## independent, as in those papers.
##
## Models, on the rows complete for every variable used:
##   mediator:  G ~ L + CG                         (linear; coefficient b1 for each column of L)
##   outcome:   T ~ L + G + CT                     (linear, logistic or Cox; t1 for L, t2 for G)
## Effects per one-unit increase in each column of L, holding the other columns fixed:
##   direct = t1, indirect = t2 * b1, total = t1 + t2 * b1
## on the mean-difference scale (continuous), log odds ratio (binary; rare outcome)
## or log hazard ratio (survival; rare outcome). The proportion mediated is
## indirect / total (continuous) or OR_D (OR_I - 1) / (OR_D OR_I - 1) (ratio scales).

cit.effects <- function(L, G, T, E = NULL, CT = NULL, CG = NULL,
                        outcome = c("continuous", "binary", "survival"), conf.level = 0.95) {
  outcome <- match.arg(outcome)
  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number between 0 and 1.", call. = FALSE)
  if (outcome == "survival" && is.null(E)) stop("survival outcomes need the event indicator 'E'.", call. = FALSE)
  if (outcome != "survival" && !is.null(E)) stop("'E' is used only with outcome = \"survival\".", call. = FALSE)
  d <- prep_data(L, G, T, CT, CG, E)
  if (ncol(d$G) != 1) stop("'G' must be a single mediator.", call. = FALSE)
  keep <- Reduce(`&`, d$ok)
  if (!all(keep)) d <- sub_data(d, keep) else d$n <- d$n.all
  n <- d$n; k <- ncol(d$L)
  if (outcome == "binary" && !all(d$T %in% c(0, 1))) stop("binary 'T' must be coded 0/1.", call. = FALSE)
  Lnames <- if (!is.null(colnames(d$L))) colnames(d$L) else if (k == 1) "L" else paste0("L", seq_len(k))

  ## mediator model
  Xm <- cbind(1, d$L, d$CG)
  qm <- qr(Xm, tol = 1e-7)
  if (qm$rank < ncol(Xm)) stop("the mediator model design is not of full rank.", call. = FALSE)
  bm <- qr.coef(qm, drop(d$G)); rm <- drop(d$G) - drop(Xm %*% bm)
  Vm <- sum(rm^2) / (n - ncol(Xm)) * chol2inv(qr.R(qm))[order(qm$pivot), order(qm$pivot)]
  ib <- 1 + seq_len(k)

  ## outcome model: columns L (1..k), G (k + 1), then CT
  Xo <- cbind(d$L, d$G, d$CT)
  if (outcome == "continuous") {
    X1 <- cbind(1, Xo); qo <- qr(X1, tol = 1e-7)
    if (qo$rank < ncol(X1)) stop("the outcome model design is not of full rank.", call. = FALSE)
    to <- qr.coef(qo, d$T); ro <- d$T - drop(X1 %*% to)
    Vo <- sum(ro^2) / (n - ncol(X1)) * chol2inv(qr.R(qo))[order(qo$pivot), order(qo$pivot)]
    to <- to[-1]; Vo <- Vo[-1, -1, drop = FALSE]
  } else if (outcome == "binary") {
    f <- suppressWarnings(glm.fit(cbind(1, Xo), d$T, family = binomial()))
    if (f$rank < ncol(Xo) + 1) stop("the outcome model design is not of full rank.", call. = FALSE)
    if (!f$converged) warning("logistic regression did not converge.", call. = FALSE)
    to <- f$coefficients[-1]
    Vo <- chol2inv(qr.R(f$qr))[order(f$qr$pivot), order(f$qr$pivot)][-1, -1, drop = FALSE]
  } else {
    f <- survival::coxph.fit(Xo, survival::Surv(d$T, d$E), strata = NULL, offset = NULL, init = NULL,
                             control = survival::coxph.control(), weights = NULL,
                             method = "efron", rownames = NULL)
    if (anyNA(f$coefficients)) stop("the outcome model design is not of full rank.", call. = FALSE)
    to <- f$coefficients; Vo <- f$var
  }
  rare <- if (outcome == "binary") mean(d$T) else if (outcome == "survival") mean(d$E) else NA
  if (!is.na(rare) && rare > 0.1)
    warning(sprintf("the %s decomposition assumes a rare outcome; here %.0f%% of observations have the event.",
                    if (outcome == "binary") "odds-ratio" else "hazard-ratio", 100 * rare), call. = FALSE)

  z <- qnorm(1 - (1 - conf.level) / 2)
  ratio <- outcome != "continuous"
  out <- lapply(seq_len(k), function(j) {
    t1 <- unname(to[j]); t2 <- unname(to[k + 1]); b1 <- unname(bm[ib[j]])
    ## parameters (t1, t2, b1) and their covariance (independent models)
    V <- matrix(0, 3, 3)
    V[1:2, 1:2] <- Vo[c(j, k + 1), c(j, k + 1)]; V[3, 3] <- Vm[ib[j], ib[j]]
    est <- c(total = t1 + t2 * b1, direct = t1, indirect = t2 * b1)
    grad <- rbind(total = c(1, b1, t2), direct = c(1, 0, 0), indirect = c(0, b1, t2))
    se <- sqrt(pmax(rowSums((grad %*% V) * grad), 0))
    pm_fun <- function(p) {
      tot <- p[1] + p[2] * p[3]; ind <- p[2] * p[3]
      if (ratio) exp(p[1]) * (exp(ind) - 1) / (exp(tot) - 1) else ind / tot
    }
    p0 <- c(t1, t2, b1); pm <- pm_fun(p0)
    h <- 1e-6 * pmax(abs(p0), 1)
    gpm <- vapply(1:3, function(i) { e <- replace(numeric(3), i, h[i]); (pm_fun(p0 + e) - pm_fun(p0 - e)) / (2 * h[i]) }, 0)
    se_pm <- sqrt(max(drop(t(gpm) %*% V %*% gpm), 0))
    lo <- est - z * se; hi <- est + z * se
    data.frame(L = Lnames[j],
               effect = c("total", "direct", "indirect", "prop.mediated"),
               scale = c(rep(switch(outcome, continuous = "difference", binary = "OR", survival = "HR"), 3), "proportion"),
               estimate = c(if (ratio) exp(est) else est, pm),
               lower = c(if (ratio) exp(lo) else lo, pm - z * se_pm),
               upper = c(if (ratio) exp(hi) else hi, pm + z * se_pm),
               se = c(se, se_pm),
               p = c(2 * pnorm(-abs(est / se)), 2 * pnorm(-abs(pm / se_pm))),
               n = n, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out); rownames(out) <- NULL
  out
}
