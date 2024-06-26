#' @name CITpkg
#' @useDynLib CITpkg
#' @importFrom Rcpp sourceCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

######################################################################
# Program Name: C_CIT_V14_CI.R
# Purpose: R CIT functions, some including C++ routines
# Programmer: Joshua Millstein
# Date: 05/03/23
#
# Input:
#   L: vector or nxp matrix of continuous instrumental variables
#   G: vector or nxp matrix of candidate causal mediators.
#   T: vector or nxp matrix of traits
#   C: vector or nxp matrix of traits
#   perm.index is n x n.perm matrix of random indices for the permutations, e.g., each column is a random permutation
#		of 1:n, where n is the number of samples and n.perm the number of permutations. For each permutation, each
#		column perm.index will be applied in therandomization approach for each component. Perm.index will preserve the
#		observed dependencies between tests in the permuted results allowing more accurate FDR confidence intervals to be computed.

#   trios: A matrix or dataframe of three columns. Each row represents a planned test to be conducted
#          and the number of rows is equal to the total number of tests. The first column is an
#          indicator for the column in L, the second is an indicator for the column in G, and the third
#          is an indicator for the column in T.
#
# Updates: 1) single continuous instrumental double variable, L, or 2) multiple instrumental double variables submitted by a matrix, L, of doubles, assuming that number of columns of L is equal to the number of L variables. 3) Allow permutation index to be added to allow dependencies between tests to be accounted for.

# If trios == NULL, then L is matrix of instrumental variables to be simultaneously included in the model, o.w. L is matrix where a single variable will be indicated by each row of trios.
##### Function to compute F test given continuous outcome and full vs reduced sets of covariates




# package pscl needed for zero inflation negative binomial

####### CIT w/ permutation results, continuous outcome, continuous L and G, possibly design matrix of L
## permutation p-values utilize permutation p-values for tests 1-3, and fstat from the parametric bootstrap.
## perm.imat is an input matrix that contains indices that specify each permutation, matrix dimenstion = sampleSize x n.perm. In order to estimate
## the over-dispersion parameter, which is necessary for estimating FDR confidence intervals, all tests
## used for the FDR estimate must be conducted under the same permutations. This is necessary to maintain
## the observed dependencies between tests.
##  under the null for test 4 (independence test)

## perm.index is n x n.perm matrix of random indices for the permutations, e.g., each column is a random permutation
##		of 1:n, where n is the number of samples and n.perm the number of permutations. For each permutation, each
##		column perm.index will be applied in therandomization approach for each component. Perm.index will preserve the
##		observed dependencies between tests in the permuted results allowing more accurate FDR confidence intervals to be computed.


# Recode NA's to -9999
ms_f = function(mat) {
  for (c_ in 1:ncol(mat)) {
    mat[is.na(mat[, c_]), c_] = -9999
  }
  return(mat)

}

# Causal Inference Test for a Binary Outcome in Version 1

#' Causal Inference Test for a Binary Outcome in Version 1
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a trait or clinical outcome. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is binary, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.bp.v1(L, G, T, C=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL)
#'
#' @param L Vector or nxp design matrix representing the instrumental variable(s).
#' @param G Continuous vector representing the potential causal mediator.
#' @param T Binary vector representing the clinical trait or outcome of interest.
#' @param C Vector or nxp design matrix representing adjustment covariates.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Seed for reproducible permutations.
#'
#' @details The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true. For component test 4, rather than using the semiparametric approach proposed by Millstein et al. (2009), here it is estimated completely by permutation, resulting in an exact test. If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.bp.v1 tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
#'
#' @return
#' A dataframe which includes the following columns:
#' \item{perm}{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
#' \item{p_cit}{CIT (omnibus) p-value}
#' \item{p_TassocL}{component p-value for the test of association between T and L.}
#' \item{p_TassocGgvnL}{component p-value for the test of association between T and G|L.}
#' \item{p_GassocLgvnT}{component p-value for the test of association between G and L|T.}
#' \item{p_LindTgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#' Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
#'
#' @author
#' Joshua Millstein, Mingzhi Ye
#'
#' @examples
#' # Sample Size
#' ss = 100
#'
#' # Errors
#' e1 = matrix(rnorm(ss), ncol=1)
#' e2 = matrix(rnorm(ss), ncol=1)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=1)
#' T = matrix(.3*G + e2, ncol=1)
#' T = ifelse( T > median(T), 1, 0 )
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests
#' results = cit.bp.v1(L, G, T)
#' results
#'
#' results = cit.bp.v1(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.bp.v1(L, G, T, C)
#' results
#'
#' results = cit.bp.v1(L, G, T, C, n.perm=5)
#' results
#'
#' @export
cit.bp.v1 = function(L,
                  G,
                  T,
                  C = NULL,
                  maxit = 10000,
                  n.perm = 0,
                  perm.index = NULL,
                  rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=100
  }
  if(is.null(perm.index) && n.perm!=0 ){
    set.seed(rseed)
    perm.index <- replicate(n.perm, sample(length(T)))
  }

  if (!is.null(perm.index)) {
    n.perm = ncol(perm.index)
    perm.index = as.matrix(perm.index)
    perm.index = perm.index - 1
  }

  if (is.vector(L)) {
    L = matrix(L, ncol = 1)
  } else {
    L = as.matrix(L)
  }
  if (is.vector(G)) {
    G = matrix(G, ncol = 1)
  } else {
    G = as.matrix(G)
  }
  if (is.vector(T)) {
    T = matrix(T, ncol = 1)
  } else {
    T = as.matrix(T)
  }
  if (!is.null(C)) {
    if (is.vector(C)) {
      C = matrix(C, ncol = 1)
    } else {
      C = as.matrix(C)
    }
  }

  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(C)) {
    aa = nrow(C) == nrow(T)
    if (!aa)
      stop("Error: rows of C must equal rows of T.")
  }

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(C))
    C = ms_f(C)
  ncolC = ncol(C)

  if (n.perm == 0) {
    aa = dim(G)[2] + dim(T)[2]
    if (aa != 2)
      stop("dim(G)[2] + dim(T)[2]  must equal 2")

    pval = 1.0
    pval1 = 1.0
    pval2 = 1.0
    pval3 = 1.0
    pval4 = 1.0 # output component p-values
    ntest = length(pval)
    nrow = dim(L)[1]
    ncol = dim(L)[2]

    if (is.null(C)) {
      citconlog2(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(nrow),
        as.integer(ncol),
        as.double(pval),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(rseed)
      )
      tmp=c(pval,pval1,pval2,pval3,pval4)

      startind = 5
    } else {
      citconlog2cvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(nrow),
        as.integer(ncol),
        as.integer(ncolC),
        as.double(pval),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(rseed)

      )
      tmp=c(pval,pval1,pval2,pval3,pval4)

      startind = 7
    } # End else is null C

    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    for (i in 1:5)
      rslts[1, i] = tmp[i]

  } else {
    # End if n.perm == 0

    if (is.null(perm.index)) {
      perm.index = matrix(NA, nrow = nrow(L), ncol = n.perm)
      for (j in 1:n.perm)
        perm.index[, j] = sample(1:nrow(L))
    }

    aa = dim(G)[2] + dim(T)[2]
    if (aa != 2)
      stop("dim(G)[2] + dim(T)[2]  must equal 2")

    trios = 0
    pval=rep(1.0, (n.perm +1))
    pval1=rep(1.0, (n.perm +1))
    pval2=rep(1.0, (n.perm +1))
    pval3=rep(1.0, (n.perm +1))
    pval4=rep(1.0, (n.perm +1))# output component p-values
    nrow = dim(L)[1]
    ncol = dim(L)[2]

    if (is.null(C) & is.null(rseed)) {
      # here permutations are not the same between multiple omnibus tests, so algorithm is slightly more computationally efficient.
      citconlog3p(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(nrow),
        as.integer(ncol),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(perm.index),
        as.integer(rseed)

      )

      startind = 3
    } else if (is.null(C)) {
      set.seed(rseed)
      citconlog3p(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(nrow),
        as.integer(ncol),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(perm.index),
        as.integer(rseed)

      )

      startind = 3
    } else {
      set.seed(rseed)
      citconlog3pcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(nrow),
        as.integer(ncol),
        as.integer(ncolC),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(perm.index),
        as.integer(rseed)

      )

      startind = 5
    } # End else is null covar and perm.imat

    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    rslts[, 1] = 0:n.perm
    rslts[, 3]=pval1
    rslts[, 4]=pval2
    rslts[, 5]=pval3
    rslts[, 6]=pval4
    for (i in 1:nrow(rslts))
      rslts[i, "p_cit"] = max(rslts[i, c("p_TassocL",
                                         "p_TassocGgvnL",
                                         "p_GassocLgvnT",
                                         "p_LindTgvnG")])

  } # End else perm > 0

  return(rslts)

} # End cit.bp.v1 function

# fdr function w/ overdispersion parameter
# In order to make CIs estimable, use the conservative approximation that at least 1 positve test was observed among permuted
fdr.od = function(obsp,
                   permp,
                   pnm,
                   ntests,
                   thres,
                   cl = 0.95,
                   od = NA) {
  z_ = qnorm(1 - (1 - cl) / 2)
  pcount = rep(NA, length(permp))
  for (p_ in 1:length(permp)) {
    permp[[p_]][, pnm] = ifelse(permp[[p_]][, pnm] <= thres,
                                1, 0)
    pcount[p_] = sum(permp[[p_]][, pnm], na.rm = TRUE)
  }
  p = mean(pcount, na.rm = TRUE) / ntests
  e_vr = ntests * p * (1 - p)
  o_vr = var(pcount, na.rm = TRUE)
  if (is.na(od)) {
    od = o_vr / e_vr
    if (!is.na(od))
      if (od < 1)
        od = 1
  }
  if (is.na(od))
    od = 1
  nperm = length(permp)
  mo = ntests
  ro = sum(obsp <= thres)
  vp = sum(pcount)
  vp1 = vp
  rslt = rep(NA, 4)
  if (ro > 0) {
    if (vp == 0)
      vp = 1
    mean.vp = vp / nperm
    fdr0 = mean.vp / ro
    pi0 = (mo - ro) / (mo - (vp / nperm))
    if (is.na(pi0))
      pi0 = 1
    if (pi0 < 0.5)
      pi0 = 0.5    # updated pi0 to limit its influence
    if (pi0 > 1)
      pi0 = 1
    fdr = fdr0 * pi0    # updated calculation of fdr to be robust to ro = mtests

    # variance of FDR
    mp = nperm * mo
    t1 = 1 / vp
    denom = mp - vp
    t2 = 1 / denom
    t3 = 1 / ro
    denom = ntests - ro
    if (denom < 1)
      denom = 1
    t4 = 1 /  denom
    s2fdr = (t1 + t2 + t3 + t4) * od
    ul = exp(log(fdr) + z_ * sqrt(s2fdr))
    ll = exp(log(fdr) - z_ * sqrt(s2fdr))

    rslt = c(fdr, ll, ul, pi0)
    rslt = ifelse(rslt > 1, 1, rslt)
    rslt = c(rslt, od, ro, vp1)
    names(rslt) = c("fdr", "fdr.ll", "fdr.ul", "pi.0", "od", "s.obs", "s.perm")
  }
  return(rslt)
} # End fdr.od



# function to combine q-values into an omnibus q-value that represents the intersection of alternative hypotheses and the union of null hypotheses
iuq = function(qvec) {
  qvec1 = 1 - qvec
  tmp = 1
  for (i in 1:length(qvec1))
    tmp = tmp * qvec1[i]
  qval = 1 - tmp
  return(qval)
} # End iuq

# wrapper function for fdr.od, gives q-values for input observed and permuted data
fdr.q.perm = function(obs.p,
                      perml,
                      pname,
                      ntests,
                      cl = .95,
                      od = NA) {
  # set q.value to minimum FDR for that p-value or larger p-values
  m = length(obs.p)
  new.order = order(obs.p)
  po = obs.p[new.order]
  qvals = rep(NA, m)
  for (tst in 1:m) {
    thresh = po[tst]
    thresh = ifelse(is.na(thresh), 1, thresh)
    thresh = ifelse(is.null(thresh), 1, thresh)
    if (thresh < 1) {
      qvals[tst] = fdr.od(obs.p,
                          perml,
                          pname,
                          ntests,
                          thresh,
                          cl = cl,
                          od = od)[1]
      qvals[1:tst] = ifelse(qvals[1:tst] > qvals[tst], qvals[tst], qvals[1:tst])
    } else
      qvals[tst] = 1
  } # End tst loop
  qvals1 = qvals[order(new.order)]
  return(qvals1)
} # End fdr.q.perm


# Millstein FDR (2013) parametric estimator, gives q-values
fdr.q.para = function(pvals) {
  # set q.value to minimum FDR for that p-value or larger p-values
  m = length(pvals)
  new.order = order(pvals)
  po = pvals[new.order]
  qvals = rep(NA, m)
  for (tst in 1:m) {
    thresh = po[tst]
    if (thresh > .99)
      qvals[tst] = 1
    if (thresh < 1) {
      S = sum(pvals <= thresh)
      Sp = m * thresh
      prod1 = Sp / S
      prod2 = (1 - S / m) / (1 - Sp / m)
      prod2 = ifelse(is.na(prod2), .5, prod2)
      prod2 = ifelse(prod2 < .5, .5, prod2)
      qvals[tst] = prod1 * prod2
    } # End if thresh
    qvals[1:tst] = ifelse(qvals[1:tst] > qvals[tst], qvals[tst], qvals[1:tst])
  } # End for tst
  qvals1 = qvals[order(new.order)]
  qvals1 = ifelse(qvals1 > 1, 1, qvals1)
  return(qvals1)
} # End fdrpara


# compute FDR qvalues from output of cit.bp or cit.cp, organized in a list with each element the output for a specific test
fdr.cit = function(cit.perm.list,
                   cl = .95,
                   c1 = NA) {
  pnms = c("p_TassocL",
           "p_TassocGgvnL",
           "p_GassocLgvnT",
           "p_LindTgvnG")
  nperm = nrow(cit.perm.list[[1]]) - 1
  ntest = length(cit.perm.list)
  perml = vector('list', nperm)
  obs = as.data.frame(matrix(NA, nrow = 0, ncol = ncol(cit.perm.list[[1]])))
  names(obs) = names(cit.perm.list[[1]])
  for (i in 1:ntest) {
    obs[i,] = cit.perm.list[[i]][1,]
    for (j in 1:nperm) {
      if (i == 1)
        perml[[j]] = obs[0,]
      perml[[j]][i,] = cit.perm.list[[i]][j + 1,]
    }
  }
  ## set 0 p-values to 1e-16
  for (pnm in pnms)
    obs[, pnm] = ifelse(obs[, pnm] < 1e-16, 1e-16, obs[, pnm])
  for (perm in 1:nperm) {
    for (pnm in pnms)
      perml[[perm]][, pnm] = ifelse(perml[[perm]][, pnm] < 1e-16, 1e-16, perml[[perm]][, pnm])
  }

  pnm.lst = vector('list', 4)
  pnm.lst[[1]] = c("q.TaL", "q.ll.TaL", "q.ul.TaL")
  pnm.lst[[2]] = c("q.TaGgvL", "q.ll.TaGgvL", "q.ul.TaGgvL")
  pnm.lst[[3]] = c("q.GaLgvT", "q.ll.GaLgvT", "q.ul.GaLgvT")
  pnm.lst[[4]] = c("q.LiTgvG", "q.ll.LiTgvG", "q.ul.LiTgvG")
  fdrmat = as.data.frame(matrix(NA, nrow = 0, ncol = 16))
  names(fdrmat) = c("p.cit",
                    "q.cit",
                    "q.cit.ll",
                    "q.cit.ul",
                    pnm.lst[[1]],
                    pnm.lst[[2]],
                    pnm.lst[[3]],
                    pnm.lst[[4]])

  for (tst in 1:nrow(obs)) {
    for (pind in 1:length(pnms)) {
      pname = pnms[pind]
      cutoff = obs[tst, pname]
      cutoff = ifelse(is.na(cutoff), 1, cutoff)
      cutoff = ifelse(is.null(cutoff), 1, cutoff)
      if (cutoff < 1) {
        fdrmat[tst, pnm.lst[[pind]]] = fdr.od(obs[, pname],
                                              perml,
                                              pname,
                                              nrow(obs),
                                              cutoff,
                                              cl = cl,
                                              od = c1)[1:3]
      } else
        fdrmat[tst, pnm.lst[[pind]]] = c(1, 1, 1)
    }
  }

  fdrmat[, pnms] = obs[, pnms]

  # p_TassocL
  op = order(fdrmat[, "p_TassocL"])
  for (tst in 1:nrow(fdrmat)) {
    aa = fdrmat[op[1:tst], "q.TaL"] > fdrmat[op[tst], "q.TaL"]
    fdrmat[op[1:tst], "q.TaL"] = ifelse(aa, fdrmat[op[tst], "q.TaL"], fdrmat[op[1:tst], "q.TaL"])
    fdrmat[op[1:tst], "q.ll.TaL"] = ifelse(aa, fdrmat[op[tst], "q.ll.TaL"], fdrmat[op[1:tst], "q.ll.TaL"])
    fdrmat[op[1:tst], "q.ul.TaL"] = ifelse(aa, fdrmat[op[tst], "q.ul.TaL"], fdrmat[op[1:tst], "q.ul.TaL"])
  }

  # p_TassocGgvnL
  op = order(fdrmat[, "p_TassocGgvnL"])
  for (tst in 1:nrow(fdrmat)) {
    aa = fdrmat[op[1:tst], "q.TaGgvL"] > fdrmat[op[tst], "q.TaGgvL"]
    fdrmat[op[1:tst], "q.TaGgvL"] = ifelse(aa, fdrmat[op[tst], "q.TaGgvL"], fdrmat[op[1:tst], "q.TaGgvL"])
    fdrmat[op[1:tst], "q.ll.TaGgvL"] = ifelse(aa, fdrmat[op[tst], "q.ll.TaGgvL"], fdrmat[op[1:tst], "q.ll.TaGgvL"])
    fdrmat[op[1:tst], "q.ul.TaGgvL"] = ifelse(aa, fdrmat[op[tst], "q.ul.TaGgvL"], fdrmat[op[1:tst], "q.ul.TaGgvL"])
  }

  # p_GassocLgvnT
  op = order(fdrmat[, "p_GassocLgvnT"])
  for (tst in 1:nrow(fdrmat)) {
    aa = fdrmat[op[1:tst], "q.GaLgvT"] > fdrmat[op[tst], "q.GaLgvT"]
    fdrmat[op[1:tst], "q.GaLgvT"] = ifelse(aa, fdrmat[op[tst], "q.GaLgvT"], fdrmat[op[1:tst], "q.GaLgvT"])
    fdrmat[op[1:tst], "q.ll.GaLgvT"] = ifelse(aa, fdrmat[op[tst], "q.ll.GaLgvT"], fdrmat[op[1:tst], "q.ll.GaLgvT"])
    fdrmat[op[1:tst], "q.ul.GaLgvT"] = ifelse(aa, fdrmat[op[tst], "q.ul.GaLgvT"], fdrmat[op[1:tst], "q.ul.GaLgvT"])
  }

  # p_LindTgvnG
  op = order(fdrmat[, "p_LindTgvnG"])
  for (tst in 1:nrow(fdrmat)) {
    aa = fdrmat[op[1:tst], "q.LiTgvG"] > fdrmat[op[tst], "q.LiTgvG"]
    fdrmat[op[1:tst], "q.LiTgvG"] = ifelse(aa, fdrmat[op[tst], "q.LiTgvG"], fdrmat[op[1:tst], "q.LiTgvG"])
    fdrmat[op[1:tst], "q.ll.LiTgvG"] = ifelse(aa, fdrmat[op[tst], "q.ll.LiTgvG"], fdrmat[op[1:tst], "q.ll.LiTgvG"])
    fdrmat[op[1:tst], "q.ul.LiTgvG"] = ifelse(aa, fdrmat[op[tst], "q.ul.LiTgvG"], fdrmat[op[1:tst], "q.ul.LiTgvG"])
  }

  # p.cit
  for (tst in 1:nrow(obs)) {
    fdrmat[tst, "p.cit"] = obs[tst, "p_cit"]
    fdrmat[tst, "q.cit"] = iuq(fdrmat[tst, c("q.TaL", "q.TaGgvL", "q.GaLgvT", "q.LiTgvG")])
    fdrmat[tst, "q.cit.ll"] = iuq(fdrmat[tst, c("q.ll.TaL", "q.ll.TaGgvL", "q.ll.GaLgvT", "q.ll.LiTgvG")])
    fdrmat[tst, "q.cit.ul"] = iuq(fdrmat[tst, c("q.ul.TaL", "q.ul.TaGgvL", "q.ul.GaLgvT", "q.ul.LiTgvG")])
  }


  op = order(fdrmat[, "p.cit"])
  for (tst in 1:nrow(fdrmat)) {
    aa = fdrmat[op[1:tst], "q.cit"] > fdrmat[op[tst], "q.cit"]
    fdrmat[op[1:tst], "q.cit"] = ifelse(aa, fdrmat[op[tst], "q.cit"], fdrmat[op[1:tst], "q.cit"])
    fdrmat[op[1:tst], "q.cit.ll"] = ifelse(aa, fdrmat[op[tst], "q.cit.ll"], fdrmat[op[1:tst], "q.cit.ll"])
    fdrmat[op[1:tst], "q.cit.ul"] = ifelse(aa, fdrmat[op[tst], "q.cit.ul"], fdrmat[op[1:tst], "q.cit.ul"])
  }

  return(fdrmat)
} # End fdr.cit






# function to test multivariate predictor vs multivariate outcome, conditioning on multivariate feature
linregM.nc = function(X, Y, C, ncp = 0) {
  X = as.matrix(X)
  Y = as.matrix(Y)
  C = as.matrix(C)
  fit = manova(Y ~ C + X)
  tmp = summary(fit, test = "Wilks")
  #  tmp$stats["X", "Pr(>F)"]
  stat.f = tmp$stats["X", "approx F"]
  df.num = tmp$stats["X", "num Df"]
  df.den = tmp$stats["X", "den Df"]
  pval.f = pf(stat.f,
              df.num,
              df.den,
              ncp = ncp,
              lower.tail = FALSE)
  return(pval.f)

} # End function linregM.nc

### Prediction approach to create a predictor of the instrumental variable using the mediator variables. Josh: Elastic Net????????????????
### The objective is to transform multiple partial mediators into a single mediator, reducing direct effects.
### Reducing direct effects will reduce bias in CIT.

# pred.m = function( L, M ){
#
#     #kvars = 10 # approximate number of final variables
#     #fit2 = glmnet( x=M, y=L, alpha=.5 )
#     #df = fit2$df[ fit2$df < (kvars+1) ]
#     #lambda = fit2$lambda[ length(df) ]
#     #pred = predict(fit2, s=lambda, newx=M)
#     #cf = as.matrix( coef(fit2, s=lambda) )
#
#     M = as.matrix(M)
#     ind = !is.na(L)
#     L = L[ind]
#     M = M[ ind, ]
#     fit.cv = cv.glmnet( x=M, y=L, alpha=.5 )
#     cf = as.matrix( coef(fit.cv, s="lambda.min") )
#     pred = predict(fit.cv, newx=M, s="lambda.min")
#
#     out.en = vector('list', 2)
#     out.en[[ 1 ]] = pred
#     out.en[[ 2 ]] = cf
#     return(out.en)
#
# } # End function pred.m


# CIT for binary outcome and permutation results, null is the empirical distribution for pvalues 1-3 and independence for p-value 4.
# Input:
#   L: vector or nxp matrix of continuous instrumental variables. If trios not equal to NULL then L includes a single instrumental variable for each test. If trios=NULL then L can be a vector, matrix or dataframe representing just one instrumental variable or alternatively a set of instrumental variables that jointly may be mediated by G.
#   G: vector or nxp matrix (if trios=NULL then G must be a single variable) of candidate causal mediators.
#   T: vector or nxp matrix of traits (if trios=NULL then T must be a single variable)
#   trios: A matrix or dataframe of three columns. Each row represents a planned test to be conducted
#          and the number of rows is equal to the total number of tests. The first column is an
#          indicator for the column in L, the second is an indicator for the column in G, and the third
#          is an indicator for the column in T. If trios not equal to NULL, then L, G, and T must be matrices or dataframes all of the same dimensions.

# apply .v2 modification, testing H2 using non-central chi-square

# Causal Inference Test for a Binary Outcome in Version 2

#' Causal Inference Test for a Binary Outcome in Version 2
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a trait or clinical outcome. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is binary, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.bp.v2(L, G, T, C=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL)
#'
#' @param L Vector or nxp design matrix representing the instrumental variable(s).
#' @param G Continuous vector representing the potential causal mediator.
#' @param T Binary vector representing the clinical trait or outcome of interest.
#' @param C Vector or nxp design matrix representing adjustment covariates.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Seed for reproducible permutations.
#'
#' @details The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true. For component test 4, rather than using the semiparametric approach proposed by Millstein et al. (2009), here it is estimated completely by permutation, resulting in an exact test. If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.bp tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
#'
#' @return
#' A dataframe which includes the following columns:
#' \item{perm}{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
#' \item{p_cit}{CIT (omnibus) p-value}
#' \item{p_TassocL}{component p-value for the test of association between T and L.}
#' \item{p_TassocGgvnL}{component p-value for the test of association between T and G|L.}
#' \item{p_GassocLgvnT}{component p-value for the test of association between G and L|T.}
#' \item{p_LindTgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#' Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
#'
#' @author
#' Joshua Millstein, Mingzhi Ye
#'
#' @examples
#' # Sample Size
#' ss = 100
#'
#' # Errors
#' e1 = matrix(rnorm(ss), ncol=1)
#' e2 = matrix(rnorm(ss), ncol=1)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=1)
#' T = matrix(.3*G + e2, ncol=1)
#' T = ifelse( T > median(T), 1, 0 )
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests
#' results = cit.bp.v2(L, G, T)
#' results
#'
#' results = cit.bp.v2(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.bp.v2(L, G, T, C)
#' results
#'
#' results = cit.bp.v2(L, G, T, C, n.perm=5)
#' results
#'
#' @export
cit.bp.v2 = function(L,
                     G,
                     T,
                     C = NULL,
                     maxit = 10000,
                     n.perm = 0,
                     perm.index = NULL,
                     rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=100
  }
  if(is.null(perm.index) && n.perm!=0 ){
    set.seed(rseed)
    perm.index <- replicate(n.perm, sample(length(T)))
  }

  if (!is.null(perm.index)) {
    n.perm = ncol(perm.index)
    perm.index = as.matrix(perm.index)
    perm.index = perm.index - 1
  }

  if (is.vector(L)) {
    L = matrix(L, ncol = 1)
  } else {
    L = as.matrix(L)
  }
  if (is.vector(G)) {
    G = matrix(G, ncol = 1)
  } else {
    G = as.matrix(G)
  }
  if (is.vector(T)) {
    T = matrix(T, ncol = 1)
  } else {
    T = as.matrix(T)
  }
  if (!is.null(C)) {
    if (is.vector(C)) {
      C = matrix(C, ncol = 1)
    } else {
      C = as.matrix(C)
    }
  }

  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(C)) {
    aa = nrow(C) == nrow(T)
    if (!aa)
      stop("Error: rows of C must equal rows of T.")
  }

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(C))
    C = ms_f(C)
  df.C = 0
  if (!is.null(C))
    df.C = ncol(C)
  n.L = dim(L)[1]
  df.L = dim(L)[2]
  pval=1.0
  pval1=1.0
  pval2=1.0
  pval3=1.0
  pval4=1.0
  pval3nc=1.0 # output component p-values

  if (n.perm == 0) {
    aa = dim(G)[2] + dim(T)[2]
    if (aa != 2)
      stop("dim(G)[2] + dim(T)[2]  must equal 2")

    if (is.null(C)) {
      citbin(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(n.L),
        as.integer(df.L),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )
      #tmp=c(pval1,pval2,pval3,pval4,pval3nc)


      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = n.L - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = n.L - (df.L + 2) # 2 is for G and intercept
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      if (fncp < 0)
        fncp = 0
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      pval3 = pf(G.p3,
                             df1 = df1,
                             df2 = df2,
                             fncp,
                             lower.tail = FALSE)

    } else {
      citbincvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(n.L),
        as.integer(df.L),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )



      # using pval3nc and df's, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = n.L - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = n.L - (df.L + 2) # 2 is for G and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      if (fncp < 0)
        fncp = 0
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      pval3 = pf(G.p3,
                             df1 = df1,
                             df2 = df2,
                             fncp,
                             lower.tail = FALSE)

    } # End else is null C

    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = pval3
    rslts[1, "p_LindTgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_LindTgvnG")])

  } else {
    # End if n.perm == 0

    aa = dim(G)[2] + dim(T)[2]
    if (aa != 2)
      stop("dim(G)[2] + dim(T)[2]  must equal 2")

    trios = 0
    pval=rep(1.0, (n.perm +1))
    pval1=rep(1.0, (n.perm +1))
    pval2=rep(1.0, (n.perm +1))
    pval3=rep(1.0, (n.perm +1))
    pval4=rep(1.0, (n.perm +1))
    pval3nc=rep(1.0, (n.perm +1)) # output component p-values
    n.L = dim(L)[1]
    df.L = dim(L)[2]

    if (is.null(C) & is.null(rseed)) {
      # here permutations are not the same between multiple omnibus tests, so algorithm is slightly more computationally efficient.
      citbinp(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(n.L),
        as.integer(df.L),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )



      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = n.L - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = n.L - (df.L + 2) # 2 is for df.T and intercept
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      fncp = ifelse(fncp < 0, 0, fncp)
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      pval3 = pf(G.p3,
                             df1 = df1,
                             df2 = df2,
                             fncp,
                             lower.tail = FALSE)

    } else if (is.null(C)) {
      set.seed(rseed)
      citbinp(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(n.L),
        as.integer(df.L),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )


      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = n.L - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = n.L - (df.L + 2) # 2 is for G and intercept
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      fncp = ifelse(fncp < 0, 0, fncp)
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      pval3 = pf(G.p3,
                             df1 = df1,
                             df2 = df2,
                             fncp,
                             lower.tail = FALSE)

    } else {
      set.seed(rseed)
      citbinpcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(n.L),
        as.integer(df.L),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )


      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = n.L - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = n.L - (df.L + 2) # 2 is for df.T and intercept
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      fncp = ifelse(fncp < 0, 0, fncp)
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      pval3 = pf(G.p3,
                             df1 = df1,
                             df2 = df2,
                             fncp,
                             lower.tail = FALSE)

    } # End else is null covar and perm.imat

    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +
                                                                                        1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = pval3[perm + 1]
      rslts[perm + 1, "p_LindTgvnG"] = pval4[perm + 1]
    }
  } # End else perm > 0

  return(rslts)

} # End cit.bp.v2 function




# Causal Inference Test for multiple mediators and a Binary Outcome in Version 1

#' Causal Inference Test for multiple mediators and a Binary Outcome in Version 1
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a trait or clinical outcome. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is binary, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.bp.m.v1(L, G, T, C=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL)
#'
#' @param L Vector or nxp design matrix representing the instrumental variable(s).
#' @param G Continuous vector representing the potential causal mediator.
#' @param T Binary vector representing the clinical trait or outcome of interest.
#' @param C Vector or nxp design matrix representing adjustment covariates.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Seed for reproducible permutations.
#'
#' @details The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true. For component test 4, rather than using the semiparametric approach proposed by Millstein et al. (2009), here it is estimated completely by permutation, resulting in an exact test. If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.bp tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
#'
#' @return
#' A dataframe which includes the following columns:
#' \item{perm}{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
#' \item{p_cit}{CIT (omnibus) p-value}
#' \item{p_TassocL}{component p-value for the test of association between T and L.}
#' \item{p_TassocGgvnL}{component p-value for the test of association between T and G|L.}
#' \item{p_GassocLgvnT}{component p-value for the test of association between G and L|T.}
#' \item{p_LindTgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#' Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
#'
#' @author
#' Joshua Millstein, Mingzhi Ye
#'
#' @examples
#' # Sample Size
#' ss = 100
#'
#' # Errors
#' e1 = matrix(rnorm(ss * 3), ncol=3)
#' e2 = matrix(rnorm(ss * 3), ncol=3)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=3)
#' T = matrix(.3*G + e2, ncol=3)
#' T <- rowSums(T)
#' T = ifelse( T > median(T), 1, 0 )
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests
#' results = cit.bp.m.v1(L, G, T)
#' results
#'
#' results = cit.bp.m.v1(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.bp.m.v1(L, G, T, C)
#' results
#'
#' results = cit.bp.m.v1(L, G, T, C, n.perm=5)
#' results
#'
#' @export
cit.bp.m.v1 = function(L,
                       G,
                       T,
                       C = NULL,
                       maxit = 10000,
                       n.perm = 0,
                       perm.index = NULL,
                       rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=100
  }
  if(is.null(perm.index) && n.perm!=0 ){
    set.seed(rseed)
    perm.index <- replicate(n.perm, sample(length(T)))
  }

  if (!is.null(perm.index)) {
    n.perm = ncol(perm.index)
    perm.index = as.matrix(perm.index)
    perm.index = perm.index - 1
  }

  if (is.vector(L)) {
    L = matrix(L, ncol = 1)
  } else {
    L = as.matrix(L)
  }
  if (is.vector(G)) {
    G = matrix(G, ncol = 1)
  } else {
    G = as.matrix(G)
  }
  if (is.vector(T)) {
    T = matrix(T, ncol = 1)
  } else {
    T = as.matrix(T)
  }
  if (!is.null(C)) {
    if (is.vector(C)) {
      C = matrix(C, ncol = 1)
    } else {
      C = as.matrix(C)
    }
  }

  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(C)) {
    aa = nrow(C) == nrow(T)
    if (!aa)
      stop("Error: rows of C must equal rows of T.")
  }

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(C))
    C = ms_f(C)

  colnames(T) = "T"
  colnames(G) = paste("G", 1:ncol(G), sep = "")
  colnames(L) = paste("L", 1:ncol(L), sep = "")
  if (!is.null(C))
    colnames(C) = paste("C", 1:ncol(C), sep = "")

  ## Remove missing
  #tmp = na.exclude(cbind(T, L, G, C))
  #if (is.null(C)) {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #} else {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #}
  #T = as.matrix(tmp[, "T"])
  #L = as.matrix(tmp[, colnames(L)])
  #G = as.matrix(tmp[, colnames(G)])
  #if (!is.null(C))
  #  C = as.matrix(tmp[, colnames(C)])
  #rm(tmp)

  df.C = 0
  if (!is.null(C))
    df.C = ncol(C)
  nobs = dim(T)[1]
  df.L = dim(L)[2]
  df.G = dim(G)[2]
  pval=1.0
  pval1=1.0
  pval2=1.0
  pval3=1.0
  pval4=1.0
  pval3nc=1.0 # output component p-values

  # if( n.resampl < n.perm ) n.resampl = n.perm
  if (!is.null(C)) {
    mydat = as.data.frame(cbind(L, G, T, C))
  } else
    mydat = as.data.frame(cbind(L, G, T))
  for (i in 1:ncol(mydat))
    mydat[, i] = as.numeric(mydat[, i])
  L.nms = paste("L", 1:ncol(L), sep = "")
  G.nms = paste("G", 1:ncol(G), sep = "")
  C.nms = NULL
  if (!is.null(C))
    C.nms = paste("C", 1:ncol(C), sep = "")
  names(mydat) = c(L.nms, G.nms, "T", C.nms)

  if (n.perm == 0) {
    if (is.null(C)) {
      citbinm(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test

      fncp = 0

      # p-value, p3: G ~ L|T
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, "T"], fncp)



    } else {
      citbinmcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test

      fncp = 0

      # p-value, p3: G ~ L|T
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, "T"], fncp)

    } # End else is null C

    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = p3
    rslts[1, "p_LindTgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_LindTgvnG")])

  } # End if n.perm == 0

  if (n.perm > 0) {
    if (is.null(perm.index)) {
      perm.index = matrix(NA, nrow = nrow(L), ncol = n.perm)
      for (j in 1:n.perm)
        perm.index[, j] = sample(1:nrow(L))
    }
    pval=rep(1.0, (n.perm +1))
    pval1=rep(1.0, (n.perm +1))
    pval2=rep(1.0, (n.perm +1))
    pval3=rep(1.0, (n.perm +1))
    pval4=rep(1.0, (n.perm +1))
    pval3nc=rep(1.0, (n.perm +1)) # output component p-values
    if (is.null(rseed))
      rseed = ceiling(runif(1) * 10000000)
    set.seed(rseed)

    if (is.null(C)) {
      citbinmp(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      for (j in 1:length(fncp)) {
        fncp[j] = 0
      }

      # p-value, p3: G ~ L|T
      p3 = rep(1, length(fncp))
      for (j in 1:length(fncp)) {
        ind.perm = 1:nrow(mydat)
        if (j > 1)
          ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, "T"], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3


    } else {
      citbinmpcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      for (j in 1:length(fncp)) {
        fncp[j] = 0
      }

      # p-value, p3: G ~ L|T
      p3 = rep(1, length(fncp))
      for (j in 1:length(fncp)) {
        ind.perm = 1:nrow(mydat)
        if (j > 1)
          ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, "T"], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3

    } # End else is null covar and perm.imat

    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = p3[perm + 1]
      rslts[perm + 1, "p_LindTgvnG"] = pval4[perm + 1]
    }
  } # End if perm > 0

  return(rslts)

} # End cit.bp.m.v1 function


# Causal Inference Test for multiple mediators and a Binary Outcome in Version 2

#' Causal Inference Test for multiple mediators and a Binary Outcome in Version 2
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a trait or clinical outcome. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is binary, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.bp.m.v2(L, G, T, C=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL)
#'
#' @param L Vector or nxp design matrix representing the instrumental variable(s).
#' @param G Continuous vector representing the potential causal mediator.
#' @param T Binary vector representing the clinical trait or outcome of interest.
#' @param C Vector or nxp design matrix representing adjustment covariates.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Seed for reproducible permutations.
#'
#' @details The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true. For component test 4, rather than using the semiparametric approach proposed by Millstein et al. (2009), here it is estimated completely by permutation, resulting in an exact test. If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.bp tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
#'
#' @return
#' A dataframe which includes the following columns:
#' \item{perm}{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
#' \item{p_cit}{CIT (omnibus) p-value}
#' \item{p_TassocL}{component p-value for the test of association between T and L.}
#' \item{p_TassocGgvnL}{component p-value for the test of association between T and G|L.}
#' \item{p_GassocLgvnT}{component p-value for the test of association between G and L|T.}
#' \item{p_LindTgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#' Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
#'
#' @author
#' Joshua Millstein, Mingzhi Ye
#'
#' @examples
#' # Sample Size
#' ss = 100
#'
#' # Errors
#' e1 = matrix(rnorm(ss * 3), ncol=3)
#' e2 = matrix(rnorm(ss * 3), ncol=3)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=3)
#' T = matrix(.3*G + e2, ncol=3)
#' T <- rowSums(T)
#' T = ifelse( T > median(T), 1, 0 )
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests
#' results = cit.bp.m.v2(L, G, T)
#' results
#'
#' results = cit.bp.m.v2(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.bp.m.v2(L, G, T, C)
#' results
#'
#' results = cit.bp.m.v2(L, G, T, C, n.perm=5)
#' results
#'
#' @export
cit.bp.m.v2 = function(L,
                       G,
                       T,
                       C = NULL,
                       maxit = 10000,
                       n.perm = 0,
                       perm.index = NULL,
                       rseed = NULL,
                       robust=TRUE) {
  permit = 1000
  if(is.null(rseed)){
    rseed=100
  }
  if(is.null(perm.index) && n.perm!=0 ){
    set.seed(rseed)
    perm.index <- replicate(n.perm, sample(length(T)))
  }

  if (!is.null(perm.index)) {
    n.perm = ncol(perm.index)
    perm.index = as.matrix(perm.index)
    perm.index = perm.index - 1
  }

  if (is.vector(L)) {
    L = matrix(L, ncol = 1)
  } else {
    L = as.matrix(L)
  }
  if (is.vector(G)) {
    G = matrix(G, ncol = 1)
  } else {
    G = as.matrix(G)
  }
  if (is.vector(T)) {
    T = matrix(T, ncol = 1)
  } else {
    T = as.matrix(T)
  }
  if (!is.null(C)) {
    if (is.vector(C)) {
      C = matrix(C, ncol = 1)
    } else {
      C = as.matrix(C)
    }
  }

  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(C)) {
    aa = nrow(C) == nrow(T)
    if (!aa)
      stop("Error: rows of C must equal rows of T.")
  }

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(C))
    C = ms_f(C)

  colnames(T) = "T"
  colnames(G) = paste("G", 1:ncol(G), sep = "")
  colnames(L) = paste("L", 1:ncol(L), sep = "")
  if (!is.null(C))
    colnames(C) = paste("C", 1:ncol(C), sep = "")

  ## Remove missing
  #tmp = na.exclude(cbind(T, L, G, C))
  #if (is.null(C)) {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #} else {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #}
  #T = as.matrix(tmp[, "T"])
  #L = as.matrix(tmp[, colnames(L)])
  #G = as.matrix(tmp[, colnames(G)])
  #if (!is.null(C))
  #  C = as.matrix(tmp[, colnames(C)])
  #rm(tmp)

  df.C = 0
  if (!is.null(C))
    df.C = ncol(C)
  nobs = dim(T)[1]
  df.L = dim(L)[2]
  df.G = dim(G)[2]
  pval=1.0
  pval1=1.0
  pval2=1.0
  pval3=1.0
  pval4=1.0
  pval3nc=1.0 # output component p-values

  # if( n.resampl < n.perm ) n.resampl = n.perm
  if (!is.null(C)) {
    mydat = as.data.frame(cbind(L, G, T, C))
  } else
    mydat = as.data.frame(cbind(L, G, T))
  for (i in 1:ncol(mydat))
    mydat[, i] = as.numeric(mydat[, i])
  L.nms = paste("L", 1:ncol(L), sep = "")
  G.nms = paste("G", 1:ncol(G), sep = "")
  C.nms = NULL
  if (!is.null(C))
    C.nms = paste("C", 1:ncol(C), sep = "")
  names(mydat) = c(L.nms, G.nms, "T", C.nms)

  if (n.perm == 0) {
    if (is.null(C)) {
      citbinm(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )



      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = nobs - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      if (fncp < 0)
        fncp = 0

      # p-value, p3: G ~ L|T
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, "T"], fncp)

    } else {
      citbinmcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )


      # using pval3nc and df's, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = nobs - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      if (fncp < 0)
        fncp = 0

      # p-value, p3: G ~ L|T
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, "T"], fncp)

    } # End else is null C

    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = p3
    rslts[1, "p_LindTgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_LindTgvnG")])

  } # End if n.perm == 0

  if (n.perm > 0) {
    if (is.null(perm.index)) {
      perm.index = matrix(NA, nrow = nrow(L), ncol = n.perm)
      for (j in 1:n.perm)
        perm.index[, j] = sample(1:nrow(L))
    }
    pval=rep(1.0, (n.perm +1))
    pval1=rep(1.0, (n.perm +1))
    pval2=rep(1.0, (n.perm +1))
    pval3=rep(1.0, (n.perm +1))
    pval4=rep(1.0, (n.perm +1))
    pval3nc=rep(1.0, (n.perm +1)) # output component p-values
    if (is.null(rseed))
      rseed = ceiling(runif(1) * 10000000)
    set.seed(rseed)

    if (is.null(C)) {
      citbinmp(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )

      # using pval3nc and df's, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = nobs - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      for (j in 1:length(fncp)) {
        if (fncp[j] < 0)
          fncp[j] = 0
      }

      # p-value, p3: G ~ L|T
      p3 = rep(1, length(fncp))
      for (j in 1:length(fncp)) {
        ind.perm = 1:nrow(mydat)
        if (j > 1)
          ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, "T"], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3

    } else {
      citbinmpcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )


      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = nobs - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for G and intercept
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      for (j in 1:length(fncp)) {
        if (fncp[j] < 0)
          fncp[j] = 0
      }

      # p-value, p3: G ~ L|T
      p3 = rep(1, length(fncp))
      for (j in 1:length(fncp)) {
        ind.perm = 1:nrow(mydat)
        if (j > 1)
          ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, "T"], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3

    } # End else is null covar and perm.imat

    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = pval3[perm + 1]
      rslts[perm + 1, "p_LindTgvnG"] = pval4[perm + 1]
    }
  } # End if perm > 0

  return(rslts)

} # End cit.bp.m.v2 function


# Causal Inference Test for a Binary Outcome

#' Causal Inference Test for a Binary Outcome
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a trait or clinical outcome. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is binary, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.bp(L, G, T, C=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL, robust=TRUE)
#'
#' @param L Vector or nxp design matrix representing the instrumental variable(s).
#' @param G Continuous vector representing the potential causal mediator.
#' @param T Binary vector representing the clinical trait or outcome of interest.
#' @param C Vector or nxp design matrix representing adjustment covariates.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Seed for reproducible permutations.
#' @param robust True for v2 algorithm and False for v1 algorithm
#'
#' @details The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true. For component test 4, rather than using the semiparametric approach proposed by Millstein et al. (2009), here it is estimated completely by permutation, resulting in an exact test. If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.bp tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
#'
#' @return
#' A dataframe which includes the following columns:
#' \item{perm}{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
#' \item{p_cit}{CIT (omnibus) p-value}
#' \item{p_TassocL}{component p-value for the test of association between T and L.}
#' \item{p_TassocGgvnL}{component p-value for the test of association between T and G|L.}
#' \item{p_GassocLgvnT}{component p-value for the test of association between G and L|T.}
#' \item{p_LindTgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#' Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
#'
#' @author
#' Joshua Millstein, Mingzhi Ye
#'
#' @examples
#' # Sample Size
#' ss = 100
#'
#' # Errors for single mediators
#' e1 = matrix(rnorm(ss), ncol=1)
#' e2 = matrix(rnorm(ss), ncol=1)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices for single mediators
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=1)
#' T = matrix(.3*G + e2, ncol=1)
#' T = ifelse( T > median(T), 1, 0 )
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests for single mediators and v1 algorithm
#' results = cit.bp(L, G, T, robust = FALSE)
#' results
#'
#' results = cit.bp(L, G, T, perm.index=perm.index, n.perm=5, robust = FALSE)
#' results
#'
#' results = cit.bp(L, G, T, C, robust = FALSE)
#' results
#'
#' results = cit.bp(L, G, T, C, n.perm=5, robust = FALSE)
#' results
#' # Run tests for single mediators and v2 algorithm
#' results = cit.bp(L, G, T)
#' results
#'
#' results = cit.bp(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.bp(L, G, T, C)
#' results
#'
#' results = cit.bp(L, G, T, C, n.perm=5)
#' results
#'
#' # Errors for multiple mediators.
#' e1 = matrix(rnorm(ss * 3), ncol=3)
#' e2 = matrix(rnorm(ss * 3), ncol=3)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices for multiple mediators
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=3)
#' T = matrix(.3*G + e2, ncol=3)
#' T <- rowSums(T)
#' T = ifelse( T > median(T), 1, 0 )
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests for multiple mediators and v1 algorithm
#' results = cit.bp(L, G, T, robust = FALSE)
#' results
#'
#' results = cit.bp(L, G, T, perm.index=perm.index, n.perm=5, robust = FALSE)
#' results
#'
#' results = cit.bp(L, G, T, C, robust = FALSE)
#' results
#'
#' results = cit.bp(L, G, T, C, n.perm=5, robust = FALSE)
#' results
#' # Run tests for multiple mediators and v2 algorithm
#' results = cit.bp(L, G, T)
#' results
#'
#' results = cit.bp(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.bp(L, G, T, C)
#' results
#'
#' results = cit.bp(L, G, T, C, n.perm=5)
#' results
#'
#' @export
cit.bp = function(L,
                  G,
                  T,
                  C = NULL,
                  maxit = 10000,
                  n.perm = 0,
                  perm.index = NULL,
                  rseed = NULL,
                  robust = TRUE
) {
  if(ncol(G) == 1){
    if(robust){
      return(cit.bp.v2(L, G, T, C, maxit, n.perm, perm.index, rseed))
    }
    else{
      return(cit.bp.v1(L, G, T, C, maxit, n.perm, perm.index, rseed))
    }
  }
  else{
    if(robust){
      return(cit.bp.m.v2(L, G, T, C, maxit, n.perm, perm.index, rseed))
    }
    else{
      return(cit.bp.m.v1(L, G, T, C, maxit, n.perm, perm.index, rseed))
    }
  }
} # End cit.bp function


# Causal Inference Test for a Continuous Outcome in Version 1

#' Causal Inference Test for a Continuous Outcome in Version 1
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a quantitative trait. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is continuous, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.cp.v1(L, G, T, C=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL)
#'
#' @param L Vector or nxp design matrix representing the instrumental variable(s).
#' @param G Continuous vector representing the potential causal mediator.
#' @param T Continuous vector representing the clinical trait or outcome of interest.
#' @param C Vector or nxp design matrix representing adjustment covariates.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Seed for reproducible permutations.
#'
#' @details Increasing maxit will increase the precision of the component test 4, the conditional independence test. This may be useful if a very small p-value is observed and high precision is desired, however, it will increase run time. The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true.  If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form, where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.cp.v1 tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
#'
#' @return
#' A dataframe which includes the following columns:
#' \item{perm}{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
#' \item{p_cit}{CIT (omnibus) p-value}
#' \item{p_TassocL}{component p-value for the test of association between T and L.}
#' \item{p_TassocGgvnL}{component p-value for the test of association between T and G|L.}
#' \item{p_GassocLgvnT}{component p-value for the test of association between G and L|T.}
#' \item{p_LindTgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#' Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
#'
#' @author
#' Joshua Millstein, Mingzhi Ye
#'
#' @examples
#' # Sample Size
#' ss = 100
#'
#' # Errors
#' e1 = matrix(rnorm(ss), ncol=1)
#' e2 = matrix(rnorm(ss), ncol=1)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=1)
#' T = matrix(.3*G + e2, ncol=1)
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests
#' results = cit.cp.v1(L, G, T)
#' results
#'
#' results = cit.cp.v1(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.cp.v1(L, G, T, C)
#' results
#'
#' results = cit.cp.v1(L, G, T, C, n.perm=5)
#' results
#'
#' @export
cit.cp.v1 = function(L,
                  G,
                  T,
                  C = NULL,
                  maxit = 10000,
                  n.perm = 0,
                  perm.index = NULL,
                  rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=100
  }
  if(is.null(perm.index) && n.perm!=0 ){
    set.seed(rseed)
    perm.index <- replicate(n.perm, sample(length(T)))
  }

  if (!is.null(perm.index)) {
    n.perm = ncol(perm.index)
    perm.index = as.matrix(perm.index)
    perm.index = perm.index - 1
  }

  if (is.vector(L)) {
    L = matrix(L, ncol = 1)
  } else {
    L = as.matrix(L)
  }
  if (is.vector(G)) {
    G = matrix(G, ncol = 1)
  } else {
    G = as.matrix(G)
  }
  if (is.vector(T)) {
    T = matrix(T, ncol = 1)
  } else {
    T = as.matrix(T)
  }
  if (!is.null(C)) {
    if (is.vector(C)) {
      C = matrix(C, ncol = 1)
    } else {
      C = as.matrix(C)
    }
  }

  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(C)) {
    aa = nrow(C) == nrow(T)
    if (!aa)
      stop("Error: rows of C must equal rows of T.")
  }

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(C))
    C = ms_f(C)
  ncolC = ncol(C)

  if (n.perm == 0) {
    aa = dim(G)[2] + dim(T)[2]
    if (aa != 2)
      stop("dim(G)[2] + dim(T)[2]  must equal 2")

    pval=1.0
    pval1=1.0
    pval2=1.0
    pval3=1.0
    pval4=1.0# output component p-values
    ntest = length(pval)
    nrow = dim(L)[1]
    ncol = dim(L)[2]

    if (is.null(C)) {
      citconlog2_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(nrow),
        as.integer(ncol),
        as.double(pval),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(rseed)

      )
      tmp=c(pval,pval1,pval2,pval3,pval4)

      startind = 5
    } else {
      citconlog2cvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(nrow),
        as.integer(ncol),
        as.integer(ncolC),
        as.double(pval),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(rseed)

      )
      tmp=c(pval,pval1,pval2,pval3,pval4)

      startind = 7
    } # End else is null C

    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    for (i in 1:5)
      rslts[1, i] = tmp[i]

  } else {
    # End if n.perm == 0

    if (is.null(perm.index)) {
      perm.index = matrix(NA, nrow = nrow(L), ncol = n.perm)
      for (j in 1:n.perm)
        perm.index[, j] = sample(1:nrow(L))
    }

    aa = dim(G)[2] + dim(T)[2]
    if (aa != 2)
      stop("dim(G)[2] + dim(T)[2]  must equal 2")

    trios = 0
    pval=rep(1.0, (n.perm +1))
    pval1=rep(1.0, (n.perm +1))
    pval2=rep(1.0, (n.perm +1))
    pval3=rep(1.0, (n.perm +1))
    pval4=rep(1.0, (n.perm +1))# output component p-values

    nrow = dim(L)[1]
    ncol = dim(L)[2]

    if (is.null(C) & is.null(rseed)) {
      # here permutations are not the same between multiple omnibus tests, so algorithm is slightly more computationally efficient.

      citconlog3p_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(nrow),
        as.integer(ncol),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(perm.index),
        as.integer(rseed)

      )

      startind = 3
    } else if (is.null(C)) {
      set.seed(rseed)
      citconlog3p_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(nrow),
        as.integer(ncol),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(perm.index),
        as.integer(rseed)

      )

      startind = 3
    } else {
      set.seed(rseed)
      citconlog3pcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(nrow),
        as.integer(ncol),
        as.integer(ncolC),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(perm.index),
        as.integer(rseed)

      )

      startind = 5
    } # End else is null covar and perm.imat

    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    rslts[, 1] = 0:n.perm
    rslts[, 3]=pval1
    rslts[, 4]=pval2
    rslts[, 5]=pval3
    rslts[, 6]=pval4
    for (i in 1:nrow(rslts))
      rslts[i, "p_cit"] = max(rslts[i, c("p_TassocL",
                                         "p_TassocGgvnL",
                                         "p_GassocLgvnT",
                                         "p_LindTgvnG")])

  } # End else perm > 0

  return(rslts)

} # End cit.cp.v1 function

# Causal Inference Test for a Continuous Outcome in Version 2

#' Causal Inference Test for a Continuous Outcome in Version 2
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a quantitative trait. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is continuous, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.cp.v2(L, G, T, C=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL)
#'
#' @param L Vector or nxp design matrix representing the instrumental variable(s).
#' @param G Continuous vector representing the potential causal mediator.
#' @param T Continuous vector representing the clinical trait or outcome of interest.
#' @param C Vector or nxp design matrix representing adjustment covariates.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Seed for reproducible permutations.
#'
#' @details Increasing maxit will increase the precision of the component test 4, the conditional independence test. This may be useful if a very small p-value is observed and high precision is desired, however, it will increase run time. The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true.  If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form, where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.cp tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
#'
#' @return
#' A dataframe which includes the following columns:
#' \item{perm}{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
#' \item{p_cit}{CIT (omnibus) p-value}
#' \item{p_TassocL}{component p-value for the test of association between T and L.}
#' \item{p_TassocGgvnL}{component p-value for the test of association between T and G|L.}
#' \item{p_GassocLgvnT}{component p-value for the test of association between G and L|T.}
#' \item{p_LindTgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#' Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
#'
#' @author
#' Joshua Millstein, Mingzhi Ye
#'
#' @examples
#' # Sample Size
#' ss = 100
#'
#' # Errors
#' e1 = matrix(rnorm(ss), ncol=1)
#' e2 = matrix(rnorm(ss), ncol=1)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=1)
#' T = matrix(.3*G + e2, ncol=1)
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests
#' results = cit.cp.v2(L, G, T)
#' results
#'
#' results = cit.cp.v2(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.cp.v2(L, G, T, C)
#' results
#'
#' results = cit.cp.v2(L, G, T, C, n.perm=5)
#' results
#'
#' @export
cit.cp.v2 = function(L,
                     G,
                     T,
                     C = NULL,
                     maxit = 10000,
                     n.perm = 0,
                     perm.index = NULL,
                     rseed = NULL) {
  permit = 1000

  if(is.null(rseed)){
    rseed=100
  }
  if(is.null(perm.index) && n.perm!=0 ){
    set.seed(rseed)
    perm.index <- replicate(n.perm, sample(length(T)))
  }

  if (!is.null(perm.index)) {
    n.perm = ncol(perm.index)
    perm.index = as.matrix(perm.index)
    perm.index = perm.index - 1
  }

  if (is.vector(L)) {
    L = matrix(L, ncol = 1)
  } else {
    L = as.matrix(L)
  }
  if (is.vector(G)) {
    G = matrix(G, ncol = 1)
  } else {
    G = as.matrix(G)
  }
  if (is.vector(T)) {
    T = matrix(T, ncol = 1)
  } else {
    T = as.matrix(T)
  }
  if (!is.null(C)) {
    if (is.vector(C)) {
      C = matrix(C, ncol = 1)
    } else {
      C = as.matrix(C)
    }
  }

  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(C)) {
    aa = nrow(C) == nrow(T)
    if (!aa)
      stop("Error: rows of C must equal rows of T.")
  }

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(C))
    C = ms_f(C)
  df.C = 0
  if (!is.null(C))
    df.C = ncol(C)
  n.L = dim(L)[1]
  df.L = dim(L)[2]

  pval=1.0
  pval1=1.0
  pval2=1.0
  pval3=1.0
  pval4=1.0
  pval3nc=1.0 # output component p-values

  if (n.perm == 0) {
    aa = dim(G)[2] + dim(T)[2]
    if (aa != 2)
      stop("dim(G)[2] + dim(T)[2]  must equal 2")

    if (is.null(C)) {
      citbin_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(n.L),
        as.integer(df.L),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )
      #tmp=c(pval1,pval2,pval3,pval4,pval3nc)


      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = n.L - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = n.L - (df.L + 2) # 2 is for G and intercept
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      if (fncp < 0)
        fncp = 0
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)# When pval3 and pval3nc are 0, G.p3, fncp, and G.nc are Inf, so eventually pval3 is Nan
      pval3 = pf(G.p3,
                 df1 = df1,
                 df2 = df2,
                 fncp,
                 lower.tail = FALSE)
      # print(c(pval3nc,fncp, G.nc, G.p3, pval3))

    } else {
      citbincvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(n.L),
        as.integer(df.L),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )



      # using pval3nc and df's, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = n.L - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = n.L - (df.L + 2) # 2 is for G and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      if (fncp < 0)
        fncp = 0
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      pval3 = pf(G.p3,
                 df1 = df1,
                 df2 = df2,
                 fncp,
                 lower.tail = FALSE)

    } # End else is null C

    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = pval3
    rslts[1, "p_LindTgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_LindTgvnG")])

  } else {
    # End if n.perm == 0

    aa = dim(G)[2] + dim(T)[2]
    if (aa != 2)
      stop("dim(G)[2] + dim(T)[2]  must equal 2")

    trios = 0

    pval=rep(1.0, (n.perm +1))
    pval1=rep(1.0, (n.perm +1))
    pval2=rep(1.0, (n.perm +1))
    pval3=rep(1.0, (n.perm +1))
    pval4=rep(1.0, (n.perm +1))
    pval3nc=rep(1.0, (n.perm +1)) # output component p-values

    n.L = dim(L)[1]
    df.L = dim(L)[2]

    if (is.null(C) & is.null(rseed)) {
      # here permutations are not the same between multiple omnibus tests, so algorithm is slightly more computationally efficient.
      citbinp_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(n.L),
        as.integer(df.L),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )



      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = n.L - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = n.L - (df.L + 2) # 2 is for df.T and intercept
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      fncp = ifelse(fncp < 0, 0, fncp)
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      pval3 = pf(G.p3,
                 df1 = df1,
                 df2 = df2,
                 fncp,
                 lower.tail = FALSE)

    } else if (is.null(C)) {
      set.seed(rseed)
      citbinp_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(n.L),
        as.integer(df.L),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )


      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = n.L - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = n.L - (df.L + 2) # 2 is for G and intercept
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      fncp = ifelse(fncp < 0, 0, fncp)
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      pval3 = pf(G.p3,
                 df1 = df1,
                 df2 = df2,
                 fncp,
                 lower.tail = FALSE)

    } else {
      set.seed(rseed)
      citbinpcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(n.L),
        as.integer(df.L),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )


      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = n.L - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = n.L - (df.L + 2) # 2 is for df.T and intercept
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      fncp = ifelse(fncp < 0, 0, fncp)
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      pval3 = pf(G.p3,
                 df1 = df1,
                 df2 = df2,
                 fncp,
                 lower.tail = FALSE)

    } # End else is null covar and perm.imat

    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +
                                                                1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = pval3[perm + 1]
      rslts[perm + 1, "p_LindTgvnG"] = pval4[perm + 1]
    }
  } # End else perm > 0

  return(rslts)

} # End cit.cp.v2 function


# Causal Inference Test for multiple mediators and a Continuous Outcome in Version 1

#' Causal Inference Test for multiple mediators and a Continuous Outcome in Version 1
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a quantitative trait. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is continuous, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.cp.m.v1(L, G, T, C=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL)
#'
#' @param L Vector or nxp design matrix representing the instrumental variable(s).
#' @param G Continuous vector representing the potential causal mediator.
#' @param T Continuous vector representing the clinical trait or outcome of interest.
#' @param C Vector or nxp design matrix representing adjustment covariates.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Seed for reproducible permutations.
#'
#' @details Increasing maxit will increase the precision of the component test 4, the conditional independence test. This may be useful if a very small p-value is observed and high precision is desired, however, it will increase run time. The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true.  If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form, where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.cp tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
#'
#' @return
#' A dataframe which includes the following columns:
#' \item{perm}{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
#' \item{p_cit}{CIT (omnibus) p-value}
#' \item{p_TassocL}{component p-value for the test of association between T and L.}
#' \item{p_TassocGgvnL}{component p-value for the test of association between T and G|L.}
#' \item{p_GassocLgvnT}{component p-value for the test of association between G and L|T.}
#' \item{p_LindTgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#' Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
#'
#' @author
#' Joshua Millstein, Mingzhi Ye
#'
#' @examples
#' # Sample Size
#' ss = 100
#'
#' # Errors
#' e1 = matrix(rnorm(ss * 3), ncol=3)
#' e2 = matrix(rnorm(ss * 3), ncol=3)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=3)
#' T = matrix(.3*G + e2, ncol=3)
#' T <- rowSums(T)
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests
#' results = cit.cp.m.v1(L, G, T)
#' results
#'
#' results = cit.cp.m.v1(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.cp.m.v1(L, G, T, C)
#' results
#'
#' results = cit.cp.m.v1(L, G, T, C, n.perm=5)
#' results
#'
#' @export
cit.cp.m.v1 = function(L,
                       G,
                       T,
                       C = NULL,
                       maxit = 10000,
                       n.perm = 0,
                       perm.index = NULL,
                       rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=100
  }
  if(is.null(perm.index) && n.perm!=0 ){
    set.seed(rseed)
    perm.index <- replicate(n.perm, sample(length(T)))
  }

  if (!is.null(perm.index)) {
    n.perm = ncol(perm.index)
    perm.index = as.matrix(perm.index)
    perm.index = perm.index - 1
  }

  if (is.vector(L)) {
    L = matrix(L, ncol = 1)
  } else {
    L = as.matrix(L)
  }
  if (is.vector(G)) {
    G = matrix(G, ncol = 1)
  } else {
    G = as.matrix(G)
  }
  if (is.vector(T)) {
    T = matrix(T, ncol = 1)
  } else {
    T = as.matrix(T)
  }
  if (!is.null(C)) {
    if (is.vector(C)) {
      C = matrix(C, ncol = 1)
    } else {
      C = as.matrix(C)
    }
  }

  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(C)) {
    aa = nrow(C) == nrow(T)
    if (!aa)
      stop("Error: rows of C must equal rows of T.")
  }


  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(C))
    C = ms_f(C)

  colnames(T) = "T"
  colnames(G) = paste("G", 1:ncol(G), sep = "")
  colnames(L) = paste("L", 1:ncol(L), sep = "")
  if (!is.null(C))
    colnames(C) = paste("C", 1:ncol(C), sep = "")

  ## Remove missing
  #tmp = na.exclude(cbind(T, L, G, C))
  #if (is.null(C)) {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #} else {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #}
  #T = as.matrix(tmp[, "T"])
  #L = as.matrix(tmp[, colnames(L)])
  #G = as.matrix(tmp[, colnames(G)])
  #if (!is.null(C))
  #  C = as.matrix(tmp[, colnames(C)])
  #rm(tmp)

  df.C = 0
  if (!is.null(C))
    df.C = ncol(C)
  nobs = dim(T)[1]
  df.L = dim(L)[2]
  df.G = dim(G)[2]

  pval=1.0
  pval1=1.0
  pval2=1.0
  pval3=1.0
  pval4=1.0
  pval3nc=1.0 # output component p-values

  # if( n.resampl < n.perm ) n.resampl = n.perm
  if (!is.null(C)) {
    mydat = as.data.frame(cbind(L, G, T, C))
  } else
    mydat = as.data.frame(cbind(L, G, T))
  for (i in 1:ncol(mydat))
    mydat[, i] = as.numeric(mydat[, i])
  L.nms = paste("L", 1:ncol(L), sep = "")
  G.nms = paste("G", 1:ncol(G), sep = "")
  C.nms = NULL
  if (!is.null(C))
    C.nms = paste("C", 1:ncol(C), sep = "")
  names(mydat) = c(L.nms, G.nms, "T", C.nms)

  if (n.perm == 0) {
    if (is.null(C)) {
      citbinm_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test

      fncp = 0

      # p-value, p3: G ~ L|T
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, "T"], fncp)

    } else {
      citbinmcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test

      fncp = 0

      # p-value, p3: G ~ L|T
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, "T"], fncp)

    } # End else is null C

    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = p3
    rslts[1, "p_LindTgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_LindTgvnG")])

  } # End if n.perm == 0

  if (n.perm > 0) {
    if (is.null(perm.index)) {
      perm.index = matrix(NA, nrow = nrow(L), ncol = n.perm)
      for (j in 1:n.perm)
        perm.index[, j] = sample(1:nrow(L))
    }
    pval=rep(1.0, (n.perm +1))
    pval1=rep(1.0, (n.perm +1))
    pval2=rep(1.0, (n.perm +1))
    pval3=rep(1.0, (n.perm +1))
    pval4=rep(1.0, (n.perm +1))
    pval3nc=rep(1.0, (n.perm +1)) # output component p-values
    if (is.null(rseed))
      rseed = ceiling(runif(1) * 10000000)
    set.seed(rseed)

    if (is.null(C)) {
      citbinmp_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      for (j in 1:length(fncp)) {
        fncp[j] = 0
      }

      # p-value, p3: G ~ L|T
      p3 = rep(1, length(fncp))
      for (j in 1:length(fncp)) {
        ind.perm = 1:nrow(mydat)
        if (j > 1)
          ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, "T"], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3

    } else {
      citbinmpcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      for (j in 1:length(fncp)) {
        fncp[j] = 0
      }

      # p-value, p3: G ~ L|T
      p3 = rep(1, length(fncp))
      for (j in 1:length(fncp)) {
        ind.perm = 1:nrow(mydat)
        if (j > 1)
          ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, "T"], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3

    } # End else is null covar and perm.imat

    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = p3[perm + 1]
      rslts[perm + 1, "p_LindTgvnG"] = pval4[perm + 1]
    }
  } # End if perm > 0

  return(rslts)

} # End cit.cp.m.v1 function


# Causal Inference Test for multiple mediators and a Continuous Outcome in Version 2

#' Causal Inference Test for multiple mediators and a Continuous Outcome in Version 2
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a quantitative trait. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is continuous, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.cp.m.v2(L, G, T, C=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL)
#'
#' @param L Vector or nxp design matrix representing the instrumental variable(s).
#' @param G Continuous vector representing the potential causal mediator.
#' @param T Continuous vector representing the clinical trait or outcome of interest.
#' @param C Vector or nxp design matrix representing adjustment covariates.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Seed for reproducible permutations.
#'
#' @details Increasing maxit will increase the precision of the component test 4, the conditional independence test. This may be useful if a very small p-value is observed and high precision is desired, however, it will increase run time. The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true.  If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form, where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.cp tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
#'
#' @return
#' A dataframe which includes the following columns:
#' \item{perm}{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
#' \item{p_cit}{CIT (omnibus) p-value}
#' \item{p_TassocL}{component p-value for the test of association between T and L.}
#' \item{p_TassocGgvnL}{component p-value for the test of association between T and G|L.}
#' \item{p_GassocLgvnT}{component p-value for the test of association between G and L|T.}
#' \item{p_LindTgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#' Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
#'
#' @author
#' Joshua Millstein, Mingzhi Ye
#'
#' @examples
#' # Sample Size
#' ss = 100
#'
#' # Errors
#' e1 = matrix(rnorm(ss * 3), ncol=3)
#' e2 = matrix(rnorm(ss * 3), ncol=3)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=3)
#' T = matrix(.3*G + e2, ncol=3)
#' T <- rowSums(T)
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests
#' results = cit.cp.m.v2(L, G, T)
#' results
#'
#' results = cit.cp.m.v2(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.cp.m.v2(L, G, T, C)
#' results
#'
#' results = cit.cp.m.v2(L, G, T, C, n.perm=5)
#' results
#'
#' @export
cit.cp.m.v2 = function(L,
                       G,
                       T,
                       C = NULL,
                       maxit = 10000,
                       n.perm = 0,
                       perm.index = NULL,
                       rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=100
  }
  if(is.null(perm.index) && n.perm!=0 ){
    set.seed(rseed)
    perm.index <- replicate(n.perm, sample(length(T)))
  }

  if (!is.null(perm.index)) {
    n.perm = ncol(perm.index)
    perm.index = as.matrix(perm.index)
    perm.index = perm.index - 1
  }

  if (is.vector(L)) {
    L = matrix(L, ncol = 1)
  } else {
    L = as.matrix(L)
  }
  if (is.vector(G)) {
    G = matrix(G, ncol = 1)
  } else {
    G = as.matrix(G)
  }
  if (is.vector(T)) {
    T = matrix(T, ncol = 1)
  } else {
    T = as.matrix(T)
  }
  if (!is.null(C)) {
    if (is.vector(C)) {
      C = matrix(C, ncol = 1)
    } else {
      C = as.matrix(C)
    }
  }

  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(C)) {
    aa = nrow(C) == nrow(T)
    if (!aa)
      stop("Error: rows of C must equal rows of T.")
  }


  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(C))
    C = ms_f(C)

  colnames(T) = "T"
  colnames(G) = paste("G", 1:ncol(G), sep = "")
  colnames(L) = paste("L", 1:ncol(L), sep = "")
  if (!is.null(C))
    colnames(C) = paste("C", 1:ncol(C), sep = "")

  ## Remove missing
  #tmp = na.exclude(cbind(T, L, G, C))
  #if (is.null(C)) {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #} else {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #}
  #T = as.matrix(tmp[, "T"])
  #L = as.matrix(tmp[, colnames(L)])
  #G = as.matrix(tmp[, colnames(G)])
  #if (!is.null(C))
  #  C = as.matrix(tmp[, colnames(C)])
  #rm(tmp)

  df.C = 0
  if (!is.null(C))
    df.C = ncol(C)
  nobs = dim(T)[1]
  df.L = dim(L)[2]
  df.G = dim(G)[2]

  pval=1.0
  pval1=1.0
  pval2=1.0
  pval3=1.0
  pval4=1.0
  pval3nc=1.0 # output component p-values

  # if( n.resampl < n.perm ) n.resampl = n.perm
  if (!is.null(C)) {
    mydat = as.data.frame(cbind(L, G, T, C))
  } else
    mydat = as.data.frame(cbind(L, G, T))
  for (i in 1:ncol(mydat))
    mydat[, i] = as.numeric(mydat[, i])
  L.nms = paste("L", 1:ncol(L), sep = "")
  G.nms = paste("G", 1:ncol(G), sep = "")
  C.nms = NULL
  if (!is.null(C))
    C.nms = paste("C", 1:ncol(C), sep = "")
  names(mydat) = c(L.nms, G.nms, "T", C.nms)

  if (n.perm == 0) {
    if (is.null(C)) {
      citbinm_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )



      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = nobs - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      if (fncp < 0)
        fncp = 0

      # p-value, p3: G ~ L|T
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, "T"], fncp)

    } else {
      citbinmcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(rseed)

      )


      # using pval3nc and df's, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = nobs - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      if (fncp < 0)
        fncp = 0

      # p-value, p3: G ~ L|T
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, "T"], fncp)

    } # End else is null C

    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = p3
    rslts[1, "p_LindTgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_LindTgvnG")])

  } # End if n.perm == 0

  if (n.perm > 0) {
    if (is.null(perm.index)) {
      perm.index = matrix(NA, nrow = nrow(L), ncol = n.perm)
      for (j in 1:n.perm)
        perm.index[, j] = sample(1:nrow(L))
    }
    pval=rep(1.0, (n.perm +1))
    pval1=rep(1.0, (n.perm +1))
    pval2=rep(1.0, (n.perm +1))
    pval3=rep(1.0, (n.perm +1))
    pval4=rep(1.0, (n.perm +1))
    pval3nc=rep(1.0, (n.perm +1)) # output component p-values
    if (is.null(rseed))
      rseed = ceiling(runif(1) * 10000000)
    set.seed(rseed)

    if (is.null(C)) {
      citbinmp_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )

      # using pval3nc and df's, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = nobs - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for df.T and intercept, covariates are not included in pval3 test
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      for (j in 1:length(fncp)) {
        if (fncp[j] < 0)
          fncp[j] = 0
      }

      # p-value, p3: G ~ L|T
      p3 = rep(1, length(fncp))
      for (j in 1:length(fncp)) {
        ind.perm = 1:nrow(mydat)
        if (j > 1)
          ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, "T"], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3

    } else {
      citbinmpcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(C),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.C),
        as.double(pval1),
        as.double(pval2),
        as.double(pval4),
        as.double(pval3nc),
        as.integer(perm.index),
        as.integer(rseed)

      )


      # using pval3nc and df, n.col, compute non-centrality parameter, lambda
      # transform pval3nc p-value to F-statistic w/ correct df, df.numerator = df.L, df.denominator = nobs - (df.L + df.T + 1), where df.T = 1
      df1 = df.L
      df2 = nobs - (df.L + 2) # 2 is for G and intercept
      G.nc = qf(pval3nc,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      fncp = G.nc * (df1 / df2) * (df2 - df1) - df1
      for (j in 1:length(fncp)) {
        if (fncp[j] < 0)
          fncp[j] = 0
      }

      # p-value, p3: G ~ L|T
      p3 = rep(1, length(fncp))
      for (j in 1:length(fncp)) {
        ind.perm = 1:nrow(mydat)
        if (j > 1)
          ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, "T"], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3

    } # End else is null covar and perm.imat

    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_LindTgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = pval3[perm + 1]
      rslts[perm + 1, "p_LindTgvnG"] = pval4[perm + 1]
    }
  } # End if perm > 0

  return(rslts)

} # End cit.cp.m.v2 function




# Causal Inference Test for a Continuous Outcome

#' Causal Inference Test for a Continuous Outcome
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a quantitative trait. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is continuous, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.cp(L, G, T, C=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL, robust=TRUE)
#'
#' @param L Vector or nxp design matrix representing the instrumental variable(s).
#' @param G Continuous vector representing the potential causal mediator.
#' @param T Continuous vector representing the clinical trait or outcome of interest.
#' @param C Vector or nxp design matrix representing adjustment covariates.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Seed for reproducible permutations.
#' @param robust True for v2 algorithm and False for v1 algorithm
#'
#' @details Increasing maxit will increase the precision of the component test 4, the conditional independence test. This may be useful if a very small p-value is observed and high precision is desired, however, it will increase run time. The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true.  If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form, where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.cp tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
#'
#' @return
#' A dataframe which includes the following columns:
#' \item{perm}{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
#' \item{p_cit}{CIT (omnibus) p-value}
#' \item{p_TassocL}{component p-value for the test of association between T and L.}
#' \item{p_TassocGgvnL}{component p-value for the test of association between T and G|L.}
#' \item{p_GassocLgvnT}{component p-value for the test of association between G and L|T.}
#' \item{p_LindTgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#' Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
#'
#' @author
#' Joshua Millstein, Mingzhi Ye
#'
#' @examples
#' # Sample Size
#' ss = 100
#'
#' # Errors for single mediators
#' e1 = matrix(rnorm(ss), ncol=1)
#' e2 = matrix(rnorm(ss), ncol=1)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices for single mediators
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=1)
#' T = matrix(.3*G + e2, ncol=1)
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests for single mediators and v1 algorithm
#' results = cit.cp(L, G, T, robust = FALSE)
#' results
#'
#' results = cit.cp(L, G, T, perm.index=perm.index, n.perm=5, robust = FALSE)
#' results
#'
#' results = cit.cp(L, G, T, C, robust = FALSE)
#' results
#'
#' results = cit.cp(L, G, T, C, n.perm=5, robust = FALSE)
#' results
#' # Run tests for single mediators and v2 algorithm
#' results = cit.cp(L, G, T)
#' results
#'
#' results = cit.cp(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.cp(L, G, T, C)
#' results
#'
#' results = cit.cp(L, G, T, C, n.perm=5)
#' results
#'
#' # Errors for multiple mediators.
#' e1 = matrix(rnorm(ss * 3), ncol=3)
#' e2 = matrix(rnorm(ss * 3), ncol=3)
#'
#' # Simulate genotypes, gene expression, covariates, and clinical trait matrices for multiple mediators
#' L = matrix(rbinom(ss*3,2,.5), ncol=3)
#' G = matrix(apply(.3*L, 1, sum) + e1, ncol=3)
#' T = matrix(.3*G + e2, ncol=3)
#' T <- rowSums(T)
#' C = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests for multiple mediators and v1 algorithm
#' results = cit.cp(L, G, T, robust = FALSE)
#' results
#'
#' results = cit.cp(L, G, T, perm.index=perm.index, n.perm=5, robust = FALSE)
#' results
#'
#' results = cit.cp(L, G, T, C, robust = FALSE)
#' results
#'
#' results = cit.cp(L, G, T, C, n.perm=5, robust = FALSE)
#' results
#' # Run tests for multiple mediators and v2 algorithm
#' results = cit.cp(L, G, T)
#' results
#'
#' results = cit.cp(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.cp(L, G, T, C)
#' results
#'
#' results = cit.cp(L, G, T, C, n.perm=5)
#' results
#'
#' @export
cit.cp = function(L,
                  G,
                  T,
                  C = NULL,
                  maxit = 10000,
                  n.perm = 0,
                  perm.index = NULL,
                  rseed = NULL,
                  robust = TRUE
) {
  if(ncol(G) == 1){
    if(robust){
      return(cit.cp.v2(L, G, T, C, maxit, n.perm, perm.index, rseed))
    }
    else{
      return(cit.cp.v1(L, G, T, C, maxit, n.perm, perm.index, rseed))
    }
  }
  else{
    if(robust){
      return(cit.cp.m.v2(L, G, T, C, maxit, n.perm, perm.index, rseed))
    }
    else{
      return(cit.cp.m.v1(L, G, T, C, maxit, n.perm, perm.index, rseed))
    }
  }
} # End cit.cp function
