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
#   L: Numeric or integer vector or nxp design matrix or dataframe representing the instrumental variable(s).
#   G: Numeric or integer vector, matrix, or dataframe representing the candidate causal mediator(s).
#   T: vector or nxp matrix of traits
#   CT: vector or nxp matrix of traits
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


cit.bp.v1 = function(L,
                  G,
                  T,
                  CT = NULL,
                  CG=NULL,
                  maxit = 10000,
                  n.perm = 0,
                  perm.index = NULL,
                  rseed = NULL) {

  permit = 1000
  if(is.null(rseed)){
    rseed=as.integer(Sys.time())
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
  if (!is.null(CT)) {
    if (is.vector(CT)) {
      CT = matrix(CT, ncol = 1)
    } else {
      CT = as.matrix(CT)
    }
  }
  if (!is.null(CG)) {
    if (is.vector(CG)) {
      CG = matrix(CG, ncol = 1)
    } else {
      CG = as.matrix(CG)
    }
  }
  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(CT)) {
    aa = nrow(CT) == nrow(T)
    if (!aa)
      stop("Error: rows of CT must equal rows of T.")
  }
  if (!is.null(CG)) {
    aa = nrow(CG) == nrow(T)
    if (!aa)
      stop("Error: rows of CG must equal rows of T.")
  }
  if (!is.null(perm.index) && nrow(L) != nrow(perm.index))
    stop("Error: rows of perm.index must equal rows of L.")
  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  ncolC = 0
  ncolCG = 0
  if (!is.null(CT)){
    CT = ms_f(CT)
    ncolC = ncol(CT)
  }
  if (!is.null(CG)){
    CG = ms_f(CG)
    ncolCG = ncol(CG)
  }
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

      citconlog2cvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.double(CG),
        as.integer(nrow),
        as.integer(ncol),
        as.integer(ncolC),
        as.integer(ncolCG),
        as.double(pval),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(rseed)

      )
      tmp=c(pval,pval1,pval2,pval3,pval4)


    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
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


      citconlog3pcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.double(CG),
        as.integer(nrow),
        as.integer(ncol),
        as.integer(ncolC),
        as.integer(ncolCG),
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


    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    rslts[, 1] = 0:n.perm
    rslts[, 3]=pval1
    rslts[, 4]=pval2
    rslts[, 5]=pval3
    rslts[, 6]=pval4
    for (i in 1:nrow(rslts))
      rslts[i, "p_cit"] = max(rslts[i, c("p_TassocL",
                                         "p_TassocGgvnL",
                                         "p_GassocLgvnT",
                                         "p_TindLgvnG")])

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
           "p_TindLgvnG")
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

  # p_TindLgvnG
  op = order(fdrmat[, "p_TindLgvnG"])
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



cit.bp.v2 = function(L,
                     G,
                     T,
                     CT = NULL,
                     CG=NULL,
                     maxit = 10000,
                     n.perm = 0,
                     perm.index = NULL,
                     rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=as.integer(Sys.time())
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
  if (!is.null(CT)) {
    if (is.vector(CT)) {
      CT = matrix(CT, ncol = 1)
    } else {
      CT = as.matrix(CT)
    }
  }
  if (!is.null(CG)) {
    if (is.vector(CG)) {
      CG = matrix(CG, ncol = 1)
    } else {
      CG = as.matrix(CG)
    }
  }
  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(CT)) {
    aa = nrow(CT) == nrow(T)
    if (!aa)
      stop("Error: rows of CT must equal rows of T.")
  }
  if (!is.null(CG)) {
    aa = nrow(CG) == nrow(T)
    if (!aa)
      stop("Error: rows of CG must equal rows of T.")
  }
  if (!is.null(perm.index) && nrow(L) != nrow(perm.index))
    stop("Error: rows of perm.index must equal rows of L.")

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(CT))
    CT = ms_f(CT)
  df.CT = 0
  if (!is.null(CT))
    df.CT = ncol(CT)
  if (!is.null(CG))
    CG = ms_f(CG)
  df.CG = 0
  if (!is.null(CG))
    df.CG = ncol(CG)
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


      citbincvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.double(CG),
        as.integer(maxit),
        as.integer(n.L),
        as.integer(df.L),
        as.integer(df.CT),
        as.integer(df.CG),
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
      fncp <- pmax(fncp, 0)
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)
      pval3 = pf(G.p3,
                             df1 = df1,
                             df2 = df2,
                             fncp,
                             lower.tail = FALSE)


    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = pval3
    rslts[1, "p_TindLgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_TindLgvnG")])

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
      citbinpcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.double(CG),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(n.L),
        as.integer(df.L),
        as.integer(df.CT),
        as.integer(df.CG),
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


    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +
                                                                                        1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = pval3[perm + 1]
      rslts[perm + 1, "p_TindLgvnG"] = pval4[perm + 1]
    }
  } # End else perm > 0

  return(rslts)

} # End cit.bp.v2 function





cit.bp.m.v1 = function(L,
                       G,
                       T,
                       CT = NULL,
                       CG=NULL,
                       maxit = 10000,
                       n.perm = 0,
                       perm.index = NULL,
                       rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=as.integer(Sys.time())
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
  if (!is.null(CT)) {
    if (is.vector(CT)) {
      CT = matrix(CT, ncol = 1)
    } else {
      CT = as.matrix(CT)
    }
  }
  if (!is.null(CG)) {
    if (is.vector(CG)) {
      CG = matrix(CG, ncol = 1)
    } else {
      CG = as.matrix(CG)
    }
  }
  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(CT)) {
    aa = nrow(CT) == nrow(T)
    if (!aa)
      stop("Error: rows of CT must equal rows of T.")
  }
  if (!is.null(CG)) {
    aa = nrow(CG) == nrow(T)
    if (!aa)
      stop("Error: rows of CG must equal rows of T.")
  }
  if (!is.null(perm.index) && nrow(L) != nrow(perm.index))
    stop("Error: rows of perm.index must equal rows of L.")

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(CT))
    CT = ms_f(CT)
  if (!is.null(CG))
    CG = ms_f(CG)
  colnames(T) = "T"
  colnames(G) = paste("G", 1:ncol(G), sep = "")
  colnames(L) = paste("L", 1:ncol(L), sep = "")
  if (!is.null(CT))
    colnames(CT) = paste("CT", 1:ncol(CT), sep = "")
  if (!is.null(CG))
    colnames(CG) = paste("CG", 1:ncol(CG), sep = "")
  ## Remove missing
  #tmp = na.exclude(cbind(T, L, G, CT))
  #if (is.null(CT)) {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #} else {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #}
  #T = as.matrix(tmp[, "T"])
  #L = as.matrix(tmp[, colnames(L)])
  #G = as.matrix(tmp[, colnames(G)])
  #if (!is.null(CT))
  #  CT = as.matrix(tmp[, colnames(CT)])
  #rm(tmp)

  df.CT = 0
  if (!is.null(CT))
    df.CT = ncol(CT)
  df.CG = 0
  if (!is.null(CG))
    df.CG = ncol(CG)
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
  mydat = as.data.frame(cbind(L, G, T, CT, CG))

  for (i in 1:ncol(mydat))
    mydat[, i] = as.numeric(mydat[, i])
  L.nms = paste("L", 1:ncol(L), sep = "")
  G.nms = paste("G", 1:ncol(G), sep = "")
  CT.nms = NULL
  if (!is.null(CT))
    CT.nms = paste("CT", 1:ncol(CT), sep = "")
  CG.nms = NULL
  if (!is.null(CG))
    CG.nms = paste("CG", 1:ncol(CG), sep = "")
  names(mydat) = c(L.nms, G.nms, "T", CT.nms, CG.nms)

  if (n.perm == 0) {
      citbinmcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.CT),
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
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, c("T", CG.nms)], fncp)



    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = p3
    rslts[1, "p_TindLgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_TindLgvnG")])

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

      citbinmpcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.CT),
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
          ind.perm = perm.index[, j - 1]+1
          #ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, c("T", CG.nms)], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3


    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = p3[perm + 1]
      rslts[perm + 1, "p_TindLgvnG"] = pval4[perm + 1]
    }
  } # End if perm > 0

  return(rslts)

} # End cit.bp.m.v1 function



cit.bp.m.v2 = function(L,
                       G,
                       T,
                       CT = NULL,
                       CG=NULL,
                       maxit = 10000,
                       n.perm = 0,
                       perm.index = NULL,
                       rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=as.integer(Sys.time())
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
  if (!is.null(CT)) {
    if (is.vector(CT)) {
      CT = matrix(CT, ncol = 1)
    } else {
      CT = as.matrix(CT)
    }
  }
  if (!is.null(CG)) {
    if (is.vector(CG)) {
      CG = matrix(CG, ncol = 1)
    } else {
      CG = as.matrix(CG)
    }
  }
  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(CT)) {
    aa = nrow(CT) == nrow(T)
    if (!aa)
      stop("Error: rows of CT must equal rows of T.")
  }
  if (!is.null(CG)) {
    aa = nrow(CG) == nrow(T)
    if (!aa)
      stop("Error: rows of CG must equal rows of T.")
  }
  if (!is.null(perm.index) && nrow(L) != nrow(perm.index))
    stop("Error: rows of perm.index must equal rows of L.")

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(CT))
    CT = ms_f(CT)
  if (!is.null(CG))
    CG = ms_f(CG)
  colnames(T) = "T"
  colnames(G) = paste("G", 1:ncol(G), sep = "")
  colnames(L) = paste("L", 1:ncol(L), sep = "")
  if (!is.null(CT))
    colnames(CT) = paste("CT", 1:ncol(CT), sep = "")
  if (!is.null(CG))
    colnames(CG) = paste("CG", 1:ncol(CG), sep = "")
  ## Remove missing
  #tmp = na.exclude(cbind(T, L, G, CT))
  #if (is.null(CT)) {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #} else {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #}
  #T = as.matrix(tmp[, "T"])
  #L = as.matrix(tmp[, colnames(L)])
  #G = as.matrix(tmp[, colnames(G)])
  #if (!is.null(CT))
  #  CT = as.matrix(tmp[, colnames(CT)])
  #rm(tmp)

  df.CT = 0
  if (!is.null(CT))
    df.CT = ncol(CT)
  df.CG = 0
  if (!is.null(CG))
    df.CG = ncol(CG)
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
  mydat = as.data.frame(cbind(L, G, T, CT, CG))

  for (i in 1:ncol(mydat))
    mydat[, i] = as.numeric(mydat[, i])
  L.nms = paste("L", 1:ncol(L), sep = "")
  G.nms = paste("G", 1:ncol(G), sep = "")
  CT.nms = NULL
  if (!is.null(CT))
    CT.nms = paste("CT", 1:ncol(CT), sep = "")
  CG.nms = NULL
  if (!is.null(CG))
    CG.nms = paste("CG", 1:ncol(CG), sep = "")
  names(mydat) = c(L.nms, G.nms, "T", CT.nms, CG.nms)

  if (n.perm == 0) {
      citbinmcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.CT),
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
      fncp <- pmax(fncp, 0)

      # p-value, p3: G ~ L|T
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, c("T", CG.nms)], fncp)


    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = p3
    rslts[1, "p_TindLgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_TindLgvnG")])

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

      citbinmpcvr(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.CT),
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
          ind.perm = perm.index[, j - 1]+1
          #ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, c("T", CG.nms)], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3


    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = pval3[perm + 1]
      rslts[perm + 1, "p_TindLgvnG"] = pval4[perm + 1]
    }
  } # End if perm > 0

  return(rslts)

} # End cit.bp.m.v2 function


# Causal Inference Test for a Binary Outcome

#' Causal Inference Test for a Binary Outcome
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to would become "would become "quantify uncertainty in a causal inference pertaining to a measured factor or factors or factors, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a trait or clinical outcome. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is binary, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.bp(L, G, T, CT=NULL, CG=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL, robust=TRUE)
#'
#' @param L Numeric or integer vector or nxp design matrix or dataframe representing the instrumental variable(s).
#' @param G Numeric or integer vector, matrix, or dataframe representing the candidate causal mediator(s).
#' @param T Binary integer vector, matrix, or dataframe representing the trait or outcome.
#' @param CT Numeric or integer vector, matrix, or dataframe representing adjustment covariates for T as the outcome.
#' @param CG Numeric or integer vector, matrix, or dataframe representing adjustment covariates for G as the outcome.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Random seed for reproducible results. The default (NULL) is to use a clock seed.
#' @param robust The robust version (the default) is a modification that precludes a positive test result (small p-value) for both the causal and reverse causal scenarios.
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
#' \item{p_TindLgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#'
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
#' CT = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#' CG = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
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
#' results = cit.bp(L, G, T, CT, CG, robust = FALSE)
#' results
#'
#' results = cit.bp(L, G, T, CT, CG, n.perm=5, robust = FALSE)
#' results
#' # Run tests for single mediators and v2 algorithm
#' results = cit.bp(L, G, T)
#' results
#'
#' results = cit.bp(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.bp(L, G, T, CT, CG)
#' results
#'
#' results = cit.bp(L, G, T, CT, CG, n.perm=5)
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
#' CT = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
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
#' results = cit.bp(L, G, T, CT, robust = FALSE)
#' results
#'
#' results = cit.bp(L, G, T, CT, n.perm=5, robust = FALSE)
#' results
#' # Run tests for multiple mediators and v2 algorithm
#' results = cit.bp(L, G, T)
#' results
#'
#' results = cit.bp(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.bp(L, G, T, CT)
#' results
#'
#' results = cit.bp(L, G, T, CT, n.perm=5)
#' results
#'
#' @export
cit.bp = function(L,
                  G,
                  T,
                  CT = NULL,
                  CG=NULL,
                  maxit = 10000,
                  n.perm = 0,
                  perm.index = NULL,
                  rseed = NULL,
                  robust = TRUE
) {
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
  if (!is.null(CT)) {
    if (is.vector(CT)) {
      CT = matrix(CT, ncol = 1)
    } else {
      CT = as.matrix(CT)
    }
  }
  if (!is.null(CG)) {
    if (is.vector(CG)) {
      CG = matrix(CG, ncol = 1)
    } else {
      CG = as.matrix(CG)
    }
  }
  if(ncol(G) == 1){
    if(robust){
      return(cit.bp.v2(L, G, T, CT, CG, maxit, n.perm, perm.index, rseed))
    }
    else{
      return(cit.bp.v1(L, G, T, CT, CG, maxit, n.perm, perm.index, rseed))
    }
  }
  else{
    if(robust){
      return(cit.bp.m.v2(L, G, T, CT, CG, maxit, n.perm, perm.index, rseed))
    }
    else{
      return(cit.bp.m.v1(L, G, T, CT, CG, maxit, n.perm, perm.index, rseed))
    }
  }
} # End cit.bp function



cit.cp.v1 = function(L,
                  G,
                  T,
                  CT = NULL,
                  CG=NULL,
                  maxit = 10000,
                  n.perm = 0,
                  perm.index = NULL,
                  rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=as.integer(Sys.time())
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
  if (!is.null(CT)) {
    if (is.vector(CT)) {
      CT = matrix(CT, ncol = 1)
    } else {
      CT = as.matrix(CT)
    }
  }
  if (!is.null(CG)) {
    if (is.vector(CG)) {
      CG = matrix(CG, ncol = 1)
    } else {
      CG = as.matrix(CG)
    }
  }

  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(CT)) {
    aa = nrow(CT) == nrow(T)
    if (!aa)
      stop("Error: rows of CT must equal rows of T.")
  }
  if (!is.null(CG)) {
    aa = nrow(CG) == nrow(T)
    if (!aa)
      stop("Error: rows of CG must equal rows of T.")
  }
  if (!is.null(perm.index) && nrow(L) != nrow(perm.index))
    stop("Error: rows of perm.index must equal rows of L.")
  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  ncolC = 0
  if (!is.null(CT)){
    CT = ms_f(CT)
    ncolC = ncol(CT)
  }
  ncolCG = 0
  if (!is.null(CG)){
    CG = ms_f(CG)
    ncolCG = ncol(CG)
  }
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

      citconlog2cvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.double(CG),
        as.integer(nrow),
        as.integer(ncol),
        as.integer(ncolC),
        as.integer(ncolCG),
        as.double(pval),
        as.double(pval1),
        as.double(pval2),
        as.double(pval3),
        as.double(pval4),
        as.integer(maxit),
        as.integer(rseed)

      )
      tmp=c(pval,pval1,pval2,pval3,pval4)


    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
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


      citconlog3pcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.double(CG),
        as.integer(nrow),
        as.integer(ncol),
        as.integer(ncolC),
        as.integer(ncolCG),
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


    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    rslts[, 1] = 0:n.perm
    rslts[, 3]=pval1
    rslts[, 4]=pval2
    rslts[, 5]=pval3
    rslts[, 6]=pval4
    for (i in 1:nrow(rslts))
      rslts[i, "p_cit"] = max(rslts[i, c("p_TassocL",
                                         "p_TassocGgvnL",
                                         "p_GassocLgvnT",
                                         "p_TindLgvnG")])

  } # End else perm > 0

  return(rslts)

} # End cit.cp.v1 function


cit.cp.v2 = function(L,
                     G,
                     T,
                     CT = NULL,
                     CG=NULL,
                     maxit = 10000,
                     n.perm = 0,
                     perm.index = NULL,
                     rseed = NULL) {
  permit = 1000

  if(is.null(rseed)){
    rseed=as.integer(Sys.time())
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
  if (!is.null(CT)) {
    if (is.vector(CT)) {
      CT = matrix(CT, ncol = 1)
    } else {
      CT = as.matrix(CT)
    }
  }
  if (!is.null(CG)) {
    if (is.vector(CG)) {
      CG = matrix(CG, ncol = 1)
    } else {
      CG = as.matrix(CG)
    }
  }
  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(CT)) {
    aa = nrow(CT) == nrow(T)
    if (!aa)
      stop("Error: rows of CT must equal rows of T.")
  }
  if (!is.null(CG)) {
    aa = nrow(CG) == nrow(T)
    if (!aa)
      stop("Error: rows of CG must equal rows of T.")
  }
  if (!is.null(perm.index) && nrow(L) != nrow(perm.index))
    stop("Error: rows of perm.index must equal rows of L.")

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(CT))
    CT = ms_f(CT)
  df.CT = 0
  if (!is.null(CT))
    df.CT = ncol(CT)
  if (!is.null(CG))
    CG = ms_f(CG)
  df.CG = 0
  if (!is.null(CG))
    df.CG = ncol(CG)
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

      citbincvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.double(CG),
        as.integer(maxit),
        as.integer(n.L),
        as.integer(df.L),
        as.integer(df.CT),
        as.integer(df.CG),
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
      fncp <- pmax(fncp, 0)
      G.p3 = qf(pval3,
                df1 = df1,
                df2 = df2,
                lower.tail = FALSE)# When pval3 and pval3nc are 0, G.p3, fncp, and G.nc are Inf, so eventually pval3 is Nan
      pval3 = pf(G.p3,
                 df1 = df1,
                 df2 = df2,
                 fncp,
                 lower.tail = FALSE)



    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = pval3
    rslts[1, "p_TindLgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_TindLgvnG")])

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

      citbinpcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.double(CG),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(n.L),
        as.integer(df.L),
        as.integer(df.CT),
        as.integer(df.CG),
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


    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +
                                                                1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = pval3[perm + 1]
      rslts[perm + 1, "p_TindLgvnG"] = pval4[perm + 1]
    }
  } # End else perm > 0

  return(rslts)

} # End cit.cp.v2 function



cit.cp.m.v1 = function(L,
                       G,
                       T,
                       CT = NULL,
                       CG=NULL,
                       maxit = 10000,
                       n.perm = 0,
                       perm.index = NULL,
                       rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=as.integer(Sys.time())
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
  if (!is.null(CT)) {
    if (is.vector(CT)) {
      CT = matrix(CT, ncol = 1)
    } else {
      CT = as.matrix(CT)
    }
  }
  if (!is.null(CG)) {
    if (is.vector(CG)) {
      CG = matrix(CG, ncol = 1)
    } else {
      CG = as.matrix(CG)
    }
  }
  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(CT)) {
    aa = nrow(CT) == nrow(T)
    if (!aa)
      stop("Error: rows of CT must equal rows of T.")
  }
  if (!is.null(CG)) {
    aa = nrow(CG) == nrow(T)
    if (!aa)
      stop("Error: rows of CG must equal rows of T.")
  }
  if (!is.null(perm.index) && nrow(L) != nrow(perm.index))
    stop("Error: rows of perm.index must equal rows of L.")

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(CT))
    CT = ms_f(CT)
  if (!is.null(CG))
    CG = ms_f(CG)
  colnames(T) = "T"
  colnames(G) = paste("G", 1:ncol(G), sep = "")
  colnames(L) = paste("L", 1:ncol(L), sep = "")
  if (!is.null(CT))
    colnames(CT) = paste("CT", 1:ncol(CT), sep = "")
  if (!is.null(CG))
    colnames(CG) = paste("CG", 1:ncol(CG), sep = "")
  ## Remove missing
  #tmp = na.exclude(cbind(T, L, G, CT))
  #if (is.null(CT)) {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #} else {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #}
  #T = as.matrix(tmp[, "T"])
  #L = as.matrix(tmp[, colnames(L)])
  #G = as.matrix(tmp[, colnames(G)])
  #if (!is.null(CT))
  #  CT = as.matrix(tmp[, colnames(CT)])
  #rm(tmp)

  df.CT = 0
  if (!is.null(CT))
    df.CT = ncol(CT)
  df.CG = 0
  if (!is.null(CG))
    df.CG = ncol(CG)
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
  mydat = as.data.frame(cbind(L, G, T, CT, CG))

  for (i in 1:ncol(mydat))
    mydat[, i] = as.numeric(mydat[, i])
  L.nms = paste("L", 1:ncol(L), sep = "")
  G.nms = paste("G", 1:ncol(G), sep = "")
  CT.nms = NULL
  if (!is.null(CT))
    CT.nms = paste("CT", 1:ncol(CT), sep = "")
  CG.nms = NULL
  if (!is.null(CG))
    CG.nms = paste("CG", 1:ncol(CG), sep = "")
  names(mydat) = c(L.nms, G.nms, "T", CT.nms, CG.nms)

  if (n.perm == 0) {
      citbinmcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.CT),
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
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, c("T", CG.nms)], fncp)


    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = p3
    rslts[1, "p_TindLgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_TindLgvnG")])

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

      citbinmpcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.CT),
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
          ind.perm = perm.index[, j - 1] + 1
          # ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, c("T", CG.nms)], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3

    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = p3[perm + 1]
      rslts[perm + 1, "p_TindLgvnG"] = pval4[perm + 1]
    }
  } # End if perm > 0

  return(rslts)

} # End cit.cp.m.v1 function



cit.cp.m.v2 = function(L,
                       G,
                       T,
                       CT = NULL,
                       CG=NULL,
                       maxit = 10000,
                       n.perm = 0,
                       perm.index = NULL,
                       rseed = NULL) {
  permit = 1000
  if(is.null(rseed)){
    rseed=as.integer(Sys.time())
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
  if (!is.null(CT)) {
    if (is.vector(CT)) {
      CT = matrix(CT, ncol = 1)
    } else {
      CT = as.matrix(CT)
    }
  }
  if (!is.null(CG)) {
    if (is.vector(CG)) {
      CG = matrix(CG, ncol = 1)
    } else {
      CG = as.matrix(CG)
    }
  }
  aa = nrow(L) == nrow(T)
  if (!aa)
    stop("Error: rows of L must equal rows of T.")
  aa = nrow(G) == nrow(T)
  if (!aa)
    stop("Error: rows of G must equal rows of T.")
  if (!is.null(CT)) {
    aa = nrow(CT) == nrow(T)
    if (!aa)
      stop("Error: rows of CT must equal rows of T.")
  }
  if (!is.null(CG)) {
    aa = nrow(CG) == nrow(T)
    if (!aa)
      stop("Error: rows of CG must equal rows of T.")
  }
  if (!is.null(perm.index) && nrow(L) != nrow(perm.index))
    stop("Error: rows of perm.index must equal rows of L.")

  L = ms_f(L)
  G = ms_f(G)
  T = ms_f(T)
  if (!is.null(CT))
    CT = ms_f(CT)
  if (!is.null(CG))
    CG = ms_f(CG)
  colnames(T) = "T"
  colnames(G) = paste("G", 1:ncol(G), sep = "")
  colnames(L) = paste("L", 1:ncol(L), sep = "")
  if (!is.null(CT))
    colnames(CT) = paste("CT", 1:ncol(CT), sep = "")
  if (!is.null(CG))
    colnames(CG) = paste("CG", 1:ncol(CG), sep = "")
  ## Remove missing
  #tmp = na.exclude(cbind(T, L, G, CT))
  #if (is.null(CT)) {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #} else {
  #  names(tmp) = c("T", colnames(L), colnames(G))
  #}
  #T = as.matrix(tmp[, "T"])
  #L = as.matrix(tmp[, colnames(L)])
  #G = as.matrix(tmp[, colnames(G)])
  #if (!is.null(CT))
  #  CT = as.matrix(tmp[, colnames(CT)])
  #rm(tmp)

  df.CT = 0
  if (!is.null(CT))
    df.CT = ncol(CT)
  df.CG = 0
  if (!is.null(CG))
    df.CG = ncol(CG)
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
  mydat = as.data.frame(cbind(L, G, T, CT, CG))

  for (i in 1:ncol(mydat))
    mydat[, i] = as.numeric(mydat[, i])
  L.nms = paste("L", 1:ncol(L), sep = "")
  G.nms = paste("G", 1:ncol(G), sep = "")
  CT.nms = NULL
  if (!is.null(CT))
    CT.nms = paste("CT", 1:ncol(CT), sep = "")
  CG.nms = NULL
  if (!is.null(CG))
    CG.nms = paste("CG", 1:ncol(CG), sep = "")
  names(mydat) = c(L.nms, G.nms, "T", CT.nms, CG.nms)

  if (n.perm == 0) {
      citbinmcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.integer(maxit),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.CT),
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
      fncp <- pmax(fncp, 0)

      # p-value, p3: G ~ L|T
      p3 = linregM.nc(mydat[, L.nms], mydat[, G.nms], mydat[, c("T", CG.nms)], fncp)


    ntest = 1
    rslts = as.data.frame(matrix(NA, nrow = ntest, ncol = 5))
    names(rslts) = c("p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    rslts[1, "p_TassocL"] = pval1
    rslts[1, "p_TassocGgvnL"] = pval2
    rslts[1, "p_GassocLgvnT"] = p3
    rslts[1, "p_TindLgvnG"] = pval4
    rslts[1, "p_cit"] = max(rslts[1, c("p_TassocL",
                                       "p_TassocGgvnL",
                                       "p_GassocLgvnT",
                                       "p_TindLgvnG")])

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

      citbinmpcvr_linear(
        as.double(L),
        as.double(G),
        as.double(T),
        as.double(CT),
        as.integer(maxit),
        as.integer(permit),
        as.integer(n.perm),
        as.integer(nobs),
        as.integer(df.L),
        as.integer(df.G),
        as.integer(df.CT),
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
          ind.perm = perm.index[, j - 1]+1
          #ind.perm = sample(1:nrow(mydat))
        tmpdat = mydat
        tmpdat[, L.nms] = mydat[ind.perm, L.nms]
        p3[j] = linregM.nc(tmpdat[, L.nms], tmpdat[, G.nms], tmpdat[, c("T", CG.nms)], fncp[j])
        rm(tmpdat)
      }
      pval3 = p3


    rslts = as.data.frame(matrix(NA, nrow = (n.perm + 1), ncol = 6))
    names(rslts) = c("perm",
                     "p_cit",
                     "p_TassocL",
                     "p_TassocGgvnL",
                     "p_GassocLgvnT",
                     "p_TindLgvnG")
    for (perm in 0:n.perm) {
      rslts[perm + 1, "perm"] = perm
      rslts[perm + 1, "p_cit"] = max(c(pval1[perm + 1], pval2[perm +1], pval3[perm + 1], pval4[perm + 1]))
      rslts[perm + 1, "p_TassocL"] = pval1[perm + 1]
      rslts[perm + 1, "p_TassocGgvnL"] = pval2[perm + 1]
      rslts[perm + 1, "p_GassocLgvnT"] = pval3[perm + 1]
      rslts[perm + 1, "p_TindLgvnG"] = pval4[perm + 1]
    }
  } # End if perm > 0

  return(rslts)

} # End cit.cp.m.v2 function




# Causal Inference Test for a Continuous Outcome

#' Causal Inference Test for a Continuous Outcome
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to would become "would become "quantify uncertainty in a causal inference pertaining to a measured factor or factors or factors, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a quantitative trait. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is continuous, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.cp(L, G, T, CT=NULL, CG=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL, robust=TRUE)
#'
#' @param L Numeric or integer vector or nxp design matrix or dataframe representing the instrumental variable(s).
#' @param G Numeric or integer vector, matrix, or dataframe representing the candidate causal mediator(s).
#' @param T Numeric or integer vector, matrix, or dataframe representing the trait or outcome.
#' @param CT Numeric or integer vector, matrix, or dataframe representing adjustment covariates for T as the outcome.
#' @param CG Numeric or integer vector, matrix, or dataframe representing adjustment covariates for G as the outcome.
#' @param maxit Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Random seed for reproducible results. The default (NULL) is to use a clock seed.
#' @param robust The robust version (the default) is a modification that precludes a positive test result (small p-value) for both the causal and reverse causal scenarios.
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
#' \item{p_TindLgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#'
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
#' CT = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#' CG = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
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
#' results = cit.cp(L, G, T, CT, CG, robust = FALSE)
#' results
#'
#' results = cit.cp(L, G, T, CT, CG, n.perm=5, robust = FALSE)
#' results
#' # Run tests for single mediators and v2 algorithm
#' results = cit.cp(L, G, T)
#' results
#'
#' results = cit.cp(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.cp(L, G, T, CT, CG)
#' results
#'
#' results = cit.cp(L, G, T, CT, CG, n.perm=5)
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
#' CT = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
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
#' results = cit.cp(L, G, T, CT, robust = FALSE)
#' results
#'
#' results = cit.cp(L, G, T, CT, n.perm=5, robust = FALSE)
#' results
#' # Run tests for multiple mediators and v2 algorithm
#' results = cit.cp(L, G, T)
#' results
#'
#' results = cit.cp(L, G, T, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.cp(L, G, T, CT)
#' results
#'
#' results = cit.cp(L, G, T, CT, n.perm=5)
#' results
#'
#' @export
cit.cp = function(L,
                  G,
                  T,
                  CT = NULL,
                  CG=NULL,
                  maxit = 10000,
                  n.perm = 0,
                  perm.index = NULL,
                  rseed = NULL,
                  robust = TRUE
) {
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
  if (!is.null(CT)) {
    if (is.vector(CT)) {
      CT = matrix(CT, ncol = 1)
    } else {
      CT = as.matrix(CT)
    }
  }
  if (!is.null(CG)) {
    if (is.vector(CG)) {
      CG = matrix(CG, ncol = 1)
    } else {
      CG = as.matrix(CG)
    }
  }
  if(ncol(G) == 1){
    if(robust){
      return(cit.cp.v2(L, G, T, CT, CG, maxit, n.perm, perm.index, rseed))
    }
    else{
      return(cit.cp.v1(L, G, T, CT, CG, maxit, n.perm, perm.index, rseed))
    }
  }
  else{
    if(robust){
      return(cit.cp.m.v2(L, G, T, CT, CG, maxit, n.perm, perm.index, rseed))
    }
    else{
      return(cit.cp.m.v1(L, G, T, CT, CG, maxit, n.perm, perm.index, rseed))
    }
  }
} # End cit.cp function


linreg = function( nms.full, nms.redu=NULL, nm.y, mydat ){

  mydat = na.exclude( mydat )

  vrs.2 = paste( nms.full, collapse="+" )
  formula2 = paste( nm.y, " ~ ", vrs.2, sep="")
  fit.full = glm( formula2 , data=mydat )

  if( is.null( nms.redu ) ){
    formula2 = paste( nm.y, " ~ 1 ", sep="")
    fit.redu  = glm( formula2 , data=mydat )
  } else {

    vrs.1 = paste( nms.redu, collapse="+" )
    formula1 = paste( nm.y, " ~ ", vrs.1, sep="")
    fit.redu  = glm( formula1 , data=mydat )

  } # End if null redu

  tmp = anova( fit.full, fit.redu, test="F" )
  pval.f = tmp$"Pr(>F)"[2]
  return( pval.f )

} # End function linreg

phreg = function( nms.full, nms.redu=NULL, nm.t, nm.e, mydat ){
  require(survival)
  mydat = na.exclude( mydat[,c(nms.full, nm.t, nm.e)] )
  my.surv = Surv( mydat[, nm.t ], mydat[, nm.e ] )
  rhs = paste( nms.full, collapse="+" )
  fmla = paste0("my.surv ~ ", rhs)
  fit.full = coxph(as.formula(fmla), data=mydat, method="efron", na.action=na.exclude)

  if( is.null( nms.redu ) ){
    fmla = "my.surv ~ 1"
    fit.redu = coxph(as.formula(fmla), data=mydat, method="efron", na.action=na.exclude)
  } else {
    rhs = paste( nms.redu, collapse="+" )
    fmla = paste0("my.surv ~ ", rhs)
    fit.redu = coxph(as.formula(fmla), data=mydat, method="efron", na.action=na.exclude)
  } # End if null redu

  tmp = anova( fit.full, fit.redu, test="chisq" )
  pval.chi = tmp$"Pr(>|Chi|)"[2]
  return( pval.chi )

} # End function phreg

## perm.index is n x n.perm matrix of random indices for the permutations, e.g., each column is a random permutation
##		of 1:n, where n is the number of samples and n.perm the number of permutations. For each permutation, each
##		column perm.index will be applied in therandomization approach for each component. Perm.index will preserve the
##		observed dependencies between tests in the permuted results allowing more accurate FDR confidence intervals to be computed.


# Causal Inference Test for survival outcomes

#' Causal Inference Test for survival outcomes
#'
#' This function implements a formal statistical hypothesis test, resulting in a p-value, to would become "would become "quantify uncertainty in a causal inference pertaining to a measured factor or factors or factors, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a quantitative trait. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is continuous, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
#'
#' @usage
#' cit.cox = function( L, G, T, E, CT=NULL, n.resampl=50, n.perm=0, perm.index=NULL, rseed=NULL )
#'
#' @param L Numeric or integer vector or nxp design matrix or dataframe representing the instrumental variable(s).
#' @param G Numeric or integer vector, matrix, or dataframe representing the candidate causal mediator(s).
#' @param T Continuous vector representing time values for time-to-event outcome.
#' @param E Vector of events, 0:censored, 1:event.
#' @param CT Vector or nxp design matrix representing adjustment covariates when Y is T.
#' @param n.resampl The number of instances of the test statistic for conditional independence (test 4) generated by ermutation (Millstein et al. 2009) under the null hypothesis of no mediation (independent effects of L on G and T). These data are used to estimate the parameters of the null distribution. The default is set to 50, which we have found to provide reasonable precision.
#' @param n.perm Number of permutations for each component test if greater than 0.
#' @param perm.index An n x n.perm matrix of permutation indices.
#' @param rseed Random seed for reproducible results. The default (NULL) is to use a clock seed.
#'
#' @details Increasing n.resampl will increase the precision of the component test 4, the conditional independence test. This may be useful if a very small p-value is observed and high precision is desired, cit.cp 7 however, it will increase run time. The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true. If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form, where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.cp tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
#'
#' @return
#' A dataframe which includes the following columns:
#' \item{perm}{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
#' \item{p_cit}{CIT (omnibus) p-value}
#' \item{p_TassocL}{component p-value for the test of association between T and L.}
#' \item{p_TassocGgvnL}{component p-value for the test of association between T and G|L.}
#' \item{p_TindLgvnG}{component p-value for the equivalence test of L ind T|G}
#'
#' @references
#' Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. PMID: 27153715.
#'
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
#' L = matrix(rbinom(ss * 3, 2, 0.5), ncol=3)
#' LS = matrix(rowSums(L), ncol = 1)
#' e1s = matrix(rowSums(e1), ncol = 1)
#' G = matrix(apply(0.3 * LS , 1, sum) + e1s, ncol=1)
#' E = rbinom(ss, 1, 0.5)
#' T = rexp(ss, rate = 0.1)
#' CT = matrix(matrix(rnorm(ss*2), ncol=1), ncol=2)
#'
#' n.perm = 5
#' perm.index = matrix(NA, nrow=ss, ncol=n.perm)
#' for(j in 1:ncol(perm.index)) perm.index[, j] = sample(1:ss)
#'
#' # Run tests
#' results = cit.cox(L, G, T, E)
#' results
#'
#' results = cit.cox(L, G, T, E, perm.index=perm.index, n.perm=5)
#' results
#'
#' results = cit.cox( L, G, T, E, CT)
#' results
#'
#' results = cit.cox( L, G, T, E, CT, n.perm=5)
#' results
#'
#' @export
cit.cox = function( L, G, T, E, CT=NULL, n.resampl=50, n.perm=0, perm.index=NULL, rseed=NULL ){

  if( !is.null(perm.index) ){
    n.perm = ncol(perm.index)
    perm.index = as.matrix( perm.index )
  }
  if( n.resampl < n.perm ) n.resampl = n.perm

  if(is.vector(L)) {
    L = as.data.frame( matrix( L, ncol=1) )
  } else {
    L = as.data.frame( as.matrix(L) )
  }
  if(is.vector(G)) {
    G = as.data.frame( matrix( G, ncol=1) )
  } else {
    G = as.data.frame( as.matrix(G) )
  }
  if(is.vector(T)) {
    T = as.data.frame( matrix( T, ncol=1) )
  } else {
    T = as.data.frame( as.matrix(T) )
  }
  if( !is.null(CT) ){
    if(is.vector(CT)) {
      CT = as.data.frame( matrix( CT, ncol=1) )
    } else {
      CT = as.data.frame( as.matrix(CT) )
    }
  }

  aa = nrow(L) == nrow(T)
  if( !aa ) stop( "Error: rows of L must equal rows of T." )
  aa = nrow(G) == nrow(T)
  if( !aa ) stop( "Error: rows of G must equal rows of T." )
  if( !is.null(CT) ){
    aa = nrow(CT) == nrow(T)
    if( !aa ) stop( "Error: rows of CT must equal rows of T." )
  }

  if( is.null(perm.index) ){
    perm.index = matrix(NA, nrow=nrow(L), ncol=n.resampl )
    for( j in 1:ncol(perm.index) ) perm.index[, j] = sample( 1:nrow(L) )
  }

  if( !is.null(CT) ){
    mydat = as.data.frame(cbind( L, G, T, E, CT ))
  } else mydat = as.data.frame(cbind( L, G, T, E ))

  for( i in 1:ncol(mydat) ){
    if(!is.element(colnames(mydat)[i], "E")){
      mydat[, i ] = as.numeric( mydat[, i ]  )
    } else mydat[, i ] = as.integer( mydat[, i ]  )
  }

  L.nms = paste("L", 1:ncol(L), sep="")
  CT.nms=NULL
  if( !is.null(CT) ) CT.nms = paste("CT", 1:ncol(CT), sep="")
  names(mydat) = c( L.nms,"G","T","E", CT.nms )

  pvec = rep(NA,3)

  # pval for T ~ L
  nm.t = "T"
  nm.e = "E"
  nms.full = c(L.nms, CT.nms)
  pvec[1] = phreg( nms.full, nms.redu=CT.nms, nm.t, nm.e, mydat )

  # pval for T ~ G|L
  nm.t = "T"
  nm.e = "E"
  nms.full = c("G", L.nms, CT.nms)
  nms.redu = c(L.nms, CT.nms)
  pvec[2] = phreg( nms.full, nms.redu, nm.t, nm.e, mydat )

  mydat1 = na.exclude(mydat)
  my.surv = Surv( mydat1[, nm.t ], mydat1[, nm.e ] )
  tmp = c( "my.surv ~ G", CT.nms )
  fmla = paste( tmp, collapse="+" )
  fit3 = coxph(as.formula(fmla), data=mydat1, method="efron")
  tmp = c( "my.surv ~ G", L.nms, CT.nms )
  fmla = paste( tmp, collapse="+" )
  fit5 = coxph(as.formula(fmla), data=mydat1, method="efron")
  chi.ind = anova(fit3,fit5)$"Chisq"[2]

  vrs.1 = paste( L.nms, collapse="+" )
  formula1 = paste( "G ~ ", vrs.1, sep="")
  fitG = lm( formula1, data=mydat, na.action=na.exclude)

  coef.g = rep(NA, length(L.nms) + 1)
  coef.g[ 1 ] = summary(fitG)$coefficients["(Intercept)",1]
  #for( i in 1:length(L.nms) ) coef.g[ i + 1 ] = summary(fitG)$coefficients[ L.nms[ i ],1]

  for( i in 1:length(L.nms) ) {
    tmp = try( summary(fitG)$coefficients[ L.nms[ i ],1], silent = TRUE )
    tmp = strsplit( as.character( tmp ), " ", fixed=TRUE )[[ 1 ]]
    coef.g[ i + 1 ] = ifelse( length( tmp ) == 1, as.numeric(tmp), 0 )
  } # End L.nms loop

  mydat[, "G.r"] = resid(fitG)

  chivecr = rep(NA,n.resampl)

  set.seed(rseed)

  for(rep in 1:n.resampl){

    if( rep <= n.perm ){
      nni  = perm.index[, rep ]
    } else {
      nni  = sample( 1:nrow(mydat) )
    }

    tmp = rep(0, nrow(mydat) )
    for( i in 1:length(L.nms) ) tmp = tmp + coef.g[ i + 1 ] * mydat[, L.nms[ i ] ]
    mydat[, "G.n"] = coef.g[ 1 ] + tmp + mydat[ nni, "G.r"]

    # chi-square for T ~ L|G.n
    mydat1 = na.exclude(mydat)
    my.surv = Surv( mydat1[, nm.t ], mydat1[, nm.e ] )
    tmp = c( "my.surv ~ G.n", CT.nms )
    fmla = paste( tmp, collapse="+" )
    fit_0 = coxph(as.formula(fmla), data=mydat1, method="efron")

    tmp = c( "my.surv ~ G.n", L.nms, CT.nms )
    fmla = paste( tmp, collapse="+" )
    fit_1 = coxph(as.formula(fmla), data=mydat1, method="efron")
    chivecr[ rep ] = anova(fit_0,fit_1)$"Chisq"[2]

  } # End rep loop

  #####Chi-square Method
  chivecr = chivecr[!is.na(chivecr)]
  df1 = anova(fit3,fit5)$Df[2]
  csncp = mean(chivecr,na.rm=TRUE) - df1
  if(csncp < 0) csncp = 0

  ######### Transform chi-square to normal
  npvals = pchisq(q=chivecr,df=df1,ncp=csncp,lower.tail=TRUE)
  nchivecr = qnorm(npvals)

  npchi = pchisq(q=chi.ind,df=df1,ncp=csncp,lower.tail=TRUE) #Transform observed chi-square statistic
  zchi = qnorm(npchi)
  pvec[3] = pnorm(zchi,mean=mean(nchivecr),sd=sd(nchivecr))

  pvalc = max(pvec)  ###Causal p-value

  pvals = c( pvalc, pvec )
  names(pvals) = c( "p_cit", "p_TassocL", "p_TassocGgvnL", "p_TindLgvnG")

  if( n.perm > 0 ){
    p.perm.ind = NA
    rep = n.resampl + 1

    if( rep <= n.perm ){
      nni  = perm.index[, rep ]
    } else {
      nni  = sample( 1:nrow(mydat) )
    }

    tmp = rep(0, nrow(mydat) )
    for( i in 1:length(L.nms) ) tmp = tmp + coef.g[ i + 1 ] * mydat[, L.nms[ i ] ]
    mydat[, "G.n"] = coef.g[ 1 ] + tmp + mydat[ nni, "G.r"]

    # chi-square for T ~ L|G.n
    mydat1 = na.exclude(mydat)
    my.surv = Surv( mydat1[, nm.t ], mydat1[, nm.e ] )
    tmp = c( "my.surv ~ G.n", CT.nms )
    fmla = paste( tmp, collapse="+" )
    fit_0 = coxph(as.formula(fmla), data=mydat1, method="efron")

    tmp = c( "my.surv ~ G.n", L.nms, CT.nms )
    fmla = paste( tmp, collapse="+" )
    fit_1 = coxph(as.formula(fmla), data=mydat1, method="efron")
    chivecr[ rep ] = anova(fit_0,fit_1)$"Chisq"[2]

    for( perm in 1:n.perm){

      chi.ind = chivecr[ perm ]
      chivecr.p = chivecr[ -perm ]
      csncp = mean(chivecr.p,na.rm=TRUE) - df1
      if(csncp < 0) csncp = 0

      ######### Transform chi-square to normal
      npvals = pchisq(q=chivecr.p,df=df1,ncp=csncp,lower.tail=TRUE)
      nchivecr = qnorm(npvals)

      npchi = pchisq(q=chi.ind,df=df1,ncp=csncp,lower.tail=TRUE) #Transform observed chi-square statistic
      zchi = qnorm(npchi)
      p.perm.ind[ perm ] = pnorm(zchi,mean=mean(nchivecr),sd=sd(nchivecr))
    } # End perm loop

    ########## permutation pvals for T ~ L, T ~ G|L
    # compute residuals and coefficients from fit
    p.perm.TasscL = NA
    p.perm.TasscGgvnL = NA

    nm.t = "T"
    nm.e = "E"

    nms.full.1 = c( L.nms, CT.nms)
    nms.redu.1 = CT.nms

    nms.full.2 = c("G", L.nms, CT.nms)
    nms.redu.2 = c(L.nms, CT.nms)

    for( perm in 1:n.perm){

      nni  = perm.index[, perm ]
      mydat.p = mydat

      mydat.p[ , L.nms ] = mydat[ nni , L.nms ]
      p.perm.TasscL[perm] = phreg( nms.full.1, nms.redu.1, nm.t, nm.e, mydat.p )
      mydat.p[ , L.nms ] = mydat[ , L.nms ]

      tmp.nms = nms.full.2[ !is.element( nms.full.2, nms.redu.2 ) ]
      mydat.p[ , tmp.nms ] = mydat[ nni , tmp.nms ]
      p.perm.TasscGgvnL[perm] = phreg( nms.full.2, nms.redu.2, nm.t, nm.e, mydat.p )
      mydat.p[ , tmp.nms ] = mydat[ , tmp.nms ]

    } # End perm loop

    rslts = as.data.frame( matrix(NA, ncol=(length(pvals) + 1) ) )
    names(rslts) = c( "perm", names(pvals) )
    rslts[ 1, "perm" ] = 0
    rslts[ 1, names(pvals) ] = pvals
    rslts[ 2:(n.perm + 1), "perm" ] = 1:n.perm
    rslts[ 2:(n.perm + 1), "p_TassocL" ] = p.perm.TasscL
    rslts[ 2:(n.perm + 1), "p_TassocGgvnL" ] = p.perm.TasscGgvnL
    rslts[ 2:(n.perm + 1), "p_TindLgvnG" ] = p.perm.ind
    for(i in 2:(n.perm+1)) rslts[ i, "p_cit" ] = max( rslts[ i, c( "p_TassocL", "p_TassocGgvnL", "p_TindLgvnG") ] )
    pvals = rslts

  } # End if n.perm

  return(pvals)

} # End function cit.cox
