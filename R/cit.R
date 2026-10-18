## User-facing functions. Argument order through `rseed` matches cit 2.3.2, so
## positional calls written for CRAN still work; `C` is kept as a deprecated
## name for `CT` (in 2.3.2 the covariates entered the outcome models only).

ct_alias <- function(CT, C) {
  if (is.null(C)) return(CT)
  if (!is.null(CT)) stop("supply covariates as 'CT' only ('C' is the old name).", call. = FALSE)
  warning("argument 'C' is deprecated; use 'CT' (covariates for the outcome models) ",
          "and optionally 'CG' (covariates for the mediator model).", call. = FALSE)
  C
}

check_args <- function(maxit, n.resampl, n.perm) {
  if (!is.numeric(maxit) || maxit < 1) stop("'maxit' must be a positive integer.", call. = FALSE)
  if (!is.numeric(n.resampl) || n.resampl < 3) stop("'n.resampl' must be at least 3.", call. = FALSE)
  if (!is.numeric(n.perm) || n.perm < 0) stop("'n.perm' must be 0 or a positive integer.", call. = FALSE)
}

cit.cp <- function(L, G, T, CT = NULL, n.resampl = 50, n.perm = 0, perm.index = NULL,
                   rseed = NULL, CG = NULL, robust = FALSE,
                   test4 = c("empirical", "parametric"), maxit = 10000, C = NULL) {
  CT <- ct_alias(CT, C); test4 <- match.arg(test4)
  check_args(maxit, n.resampl, n.perm)
  cit_core("linear", L, G, T, CT = CT, CG = CG, robust = robust, test4 = test4,
           maxit = maxit, n.resampl = n.resampl, n.perm = n.perm,
           perm.index = perm.index, rseed = rseed)
}

cit.bp <- function(L, G, T, CT = NULL, maxit = 10000, n.perm = 0, perm.index = NULL,
                   rseed = NULL, CG = NULL, robust = FALSE,
                   test4 = c("empirical", "parametric"), n.resampl = 50, C = NULL) {
  CT <- ct_alias(CT, C); test4 <- match.arg(test4)
  check_args(maxit, n.resampl, n.perm)
  res <- cit_core("logistic", L, G, T, CT = CT, CG = CG, robust = robust, test4 = test4,
                  maxit = maxit, n.resampl = n.resampl, n.perm = n.perm,
                  perm.index = perm.index, rseed = rseed)
  # as in 2.3.2, cit.bp returns a one-row data frame when n.perm = 0
  if (is.data.frame(res)) res else as.data.frame(as.list(res))
}

cit.cox <- function(L, G, T, E, CT = NULL, n.resampl = 50, n.perm = 0, perm.index = NULL,
                    rseed = NULL, CG = NULL, test4 = c("empirical", "parametric"),
                    maxit = 10000) {
  test4 <- match.arg(test4)
  check_args(maxit, n.resampl, n.perm)
  cit_core("cox", L, G, T, E = E, CT = CT, CG = CG, robust = FALSE, test4 = test4,
           maxit = maxit, n.resampl = n.resampl, n.perm = n.perm,
           perm.index = perm.index, rseed = rseed)
}
