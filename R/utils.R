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

## Restrict permutations of all n rows to the complete-case rows `keep`.
## Where a kept row's source row is also kept, the pairing is preserved; the
## remaining kept rows receive the unused kept sources in order. With no
## missing data this returns perm.index unchanged, so permutations stay
## aligned across tests that share a perm.index.
restrict_perm <- function(perm.index, keep) {
  idx <- which(keep)
  apply(perm.index, 2, function(p) {
    src <- p[idx]
    ok <- keep[src]
    src[!ok] <- setdiff(p[keep[p]], src[ok])
    match(src, idx)
  })
}

## Assemble and clean the inputs: matrices, row checks, complete cases.
prep_data <- function(L, G, T, CT, CG, E = NULL) {
  L <- as_dmat(L, "L"); G <- as_dmat(G, "G"); T <- as_dmat(T, "T")
  CT <- as_dmat(CT, "CT"); CG <- as_dmat(CG, "CG"); E <- as_dmat(E, "E")
  n <- nrow(T)
  if (ncol(T) != 1) stop("'T' must be a single variable.", call. = FALSE)
  for (nm in c("L", "G", "CT", "CG", "E")) {
    x <- get(nm)
    if (!is.null(x) && nrow(x) != n) stop(sprintf("rows of %s must equal rows of T.", nm), call. = FALSE)
  }
  keep <- stats::complete.cases(L, G, T, CT, CG, E)
  sub <- function(x) if (is.null(x)) NULL else x[keep, , drop = FALSE]
  list(L = sub(L), G = sub(G), T = drop(sub(T)), CT = sub(CT), CG = sub(CG),
       E = if (is.null(E)) NULL else drop(sub(E)), keep = keep, n.all = n, n = sum(keep))
}

## columns of X whose pivot falls within the numerical rank of cbind(base, X)
independent_cols <- function(base, X) {
  q <- qr(cbind(base, X), tol = 1e-7)
  piv <- q$pivot[seq_len(q$rank)] - ncol(base)
  sort(piv[piv > 0])
}
