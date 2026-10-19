## Input handling shared by cit.cp, cit.bp and cit.cox.

as_dmat <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) x <- data.matrix(x)
  x <- as.matrix(x)
  if (!is.numeric(x) && !is.logical(x)) stop(sprintf("'%s' must be numeric.", name), call. = FALSE)
  storage.mode(x) <- "double"
  x
}

## Reproducible randomness: with rseed, use that seed and restore the caller's
## RNG state afterwards; without it, draw from the caller's stream so that
## set.seed() before the call controls the result.
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = genv)
          else if (exists(".Random.seed", envir = genv, inherits = FALSE))
            rm(".Random.seed", envir = genv))
  set.seed(seed)
  expr
}

## Check a user-supplied perm.index (n rows, each column a permutation of 1:n)
check_perm_index <- function(perm.index, n) {
  perm.index <- as.matrix(perm.index)
  if (nrow(perm.index) != n) stop("'perm.index' must have one row per observation.", call. = FALSE)
  ok <- apply(perm.index, 2, function(p) !anyNA(p) && all(sort(p) == seq_len(n)))
  if (!all(ok)) stop("each column of 'perm.index' must be a permutation of 1:n.", call. = FALSE)
  storage.mode(perm.index) <- "integer"
  perm.index
}

## Restrict permutations of all n rows to the rows `keep` (used for test 4,
## whose resampled residuals exist only for its complete rows). Where a kept
## row's source row is also kept, the pairing is preserved; the remaining kept
## rows receive the unused kept sources in order. With no missing data this
## returns perm.index unchanged.
restrict_perm <- function(perm.index, keep) {
  idx <- which(keep)
  apply(perm.index, 2, function(p) {
    src <- p[idx]
    ok <- keep[src]
    src[!ok] <- setdiff(p[keep[p]], src[ok])
    match(src, idx)
  })
}

## Assemble and check the inputs. Rows are not dropped here: each component
## test uses the rows complete for its own variables (see run_cit), so `ok`
## records which rows are complete for each variable.
prep_data <- function(L, G, T, CT, CG, E = NULL) {
  L <- as_dmat(L, "L"); G <- as_dmat(G, "G"); T <- as_dmat(T, "T")
  CT <- as_dmat(CT, "CT"); CG <- as_dmat(CG, "CG"); E <- as_dmat(E, "E")
  n <- nrow(T)
  if (ncol(T) != 1) stop("'T' must be a single variable.", call. = FALSE)
  for (nm in c("L", "G", "CT", "CG", "E")) {
    x <- get(nm)
    if (!is.null(x) && nrow(x) != n) stop(sprintf("rows of %s must equal rows of T.", nm), call. = FALSE)
  }
  cc <- function(x) if (is.null(x)) rep(TRUE, n) else stats::complete.cases(x)
  ok <- list(L = cc(L), G = cc(G), T = cc(T) & cc(E), CT = cc(CT), CG = cc(CG))
  list(L = L, G = G, T = drop(T), CT = CT, CG = CG, E = if (is.null(E)) NULL else drop(E),
       ok = ok, n.all = n)
}

## the data restricted to rows `r`, with instrument Lx (all n rows)
sub_data <- function(d, r, Lx = d$L) {
  sub <- function(x) if (is.null(x)) NULL else x[r, , drop = FALSE]
  list(L = sub(Lx), G = sub(d$G), T = d$T[r], CT = sub(d$CT), CG = sub(d$CG),
       E = if (is.null(d$E)) NULL else d$E[r], n = sum(r))
}

## columns of X whose pivot falls within the numerical rank of cbind(base, X)
independent_cols <- function(base, X) {
  q <- qr(cbind(base, X), tol = 1e-7)
  piv <- q$pivot[seq_len(q$rank)] - ncol(base)
  sort(piv[piv > 0])
}
