\name{fdr.cit}
\alias{fdr.cit}
\title{
  Omnibus FDR Values for CIT
}
\description{
   Results including permutation results from \code{\link{cit.cp}}, \code{\link{cit.bp}} or \code{\link{cit.cox}} from multiple tests are used in a permutation-based approach to compute false discovery rates (q-values) for each test. q-values are computed for the omnibus test as well as all component tests. q-values are returned with confidence intervals to quantify uncertainty in the estimate. 
}
\usage{
   fdr.cit( cit.perm.list, cl=.95, c1=NA )
}

\arguments{
  \item{cit.perm.list}{
     List, where each item is a dataframe output from cit.cp, cit.bp or cit.cox, which must be run with n.perm > 0 and the same perm.index. 
}
  \item{cl}{
     Confidence level for the q-value confidence interval (default is .95).
}
  \item{c1}{
     Over-dispersion parameter. Setting the over-dispersion parameter to one corresponds to the assumption that all tests (referring to multiple omnibus tests conducted using different sets of variables) are independent. The default (c1=NA) is to use an empirically estimated over-dispersion parameter based on the permutation results.
}
}
\details{
  For the initial cit.cp and cit.bp test computations, we suggest setting n.perm to 100 or greater for moderate numbers of tests, say 10s or 100s. Larger numbers of permutations will generate more precise FDR estimates but require substantial computation time. If a large number of tests, say 1000s, have been conducted, then n.perm can be as small as 10 and still generate good FDR estimates and confidence intervals. FDR confidence intervals are especially useful to quantify uncertainty when FDR value > 0.05. The over-dispersion parameter should be set to one if it is known that all tests are independent, however, this is rarely the case in typical multiple testing settings. The fdr.cit function can also be used with simulated data in order to generate power calculations. In this case the over-dispersion parameter may be estimated from pilot data or other datasets and fixed to the prior estimate using the c1 argument. The omnibus CIT FDR value combines the component test FDR values (four component tests, or three for cit.cox with test3 = FALSE) in an intersection-union type of approach analogous to the method used to generate the omnibus CIT p-value; see below.
  
  \emph{The omnibus q-value.} The omnibus q-value is a q-value for the omnibus p-value \eqn{p_{cit} = \max_k p_k}{p_cit = max_k p_k}, built from the component tests. For a test with omnibus p-value \eqn{t}, each component FDR is estimated by permutation (Millstein and Volfson 2013) at the common threshold \eqn{t}, \eqn{FDR_k(t)}, and the estimates are combined as \eqn{Q(t) = 1 - \prod_k (1 - FDR_k(t))}{Q(t) = 1 - prod_k (1 - FDR_k(t))}; q.cit is the minimum of \eqn{Q} over the omnibus p-values \eqn{\ge t}{>= t}. The rejection regions \eqn{\{p_{cit} \le t\}}{\{p_cit <= t\}} are nested, so selecting the tests with q.cit \eqn{\le \alpha}{<= alpha} selects the tests with \eqn{p_{cit}}{p_cit} below a single cutoff, as with ordinary q-values (Storey 2002, 2003). (Versions before 3.0.0 evaluated each component at the test's own component p-value. Those regions differ between tests and are not nested, and selecting q.cit \eqn{\le \alpha}{<= alpha} did not in general control the FDR. The omnibus q-values are never smaller than in 2.3.2, up to permutation step effects. The component q-values are unchanged.)

  The combination is justified under the two-groups model (Efron et al. 2001; Efron 2008), in which each hypothesis is null or non-null at random and its p-value is drawn given that status; under this model the positive FDR (pFDR) of a rejection region equals the probability that a hypothesis in the region is null (Storey 2003). Let \eqn{H = (H_1, \ldots, H_K)}{H = (H_1, ..., H_K)} be the null statuses of the K component hypotheses of a test. A CIT call is a false discovery when at least one component null holds (the composite null). Suppose (1) over the tests analysed together, the K component null statuses are mutually independent, and (2) given \eqn{H}, the component p-values are independent, with each \eqn{p_k} depending only on its own status \eqn{H_k}. Then \eqn{1 - \prod_k (1 - pFDR_k(t))}{1 - prod_k (1 - pFDR_k(t))} equals the pFDR of the call set \eqn{\{p_{cit} \le t\}}{\{p_cit <= t\}}, and q.cit estimates the corresponding q-value consistently when the component estimates are consistent. These conditions concern how component null statuses co-occur among the tests analysed together, which depends on how those tests were selected, and they can fail. If many tests are null for one component and strongly non-null for the others, the omnibus q-value can understate the FDR; if tests under the composite null tend to be null for every component, it can overstate it.

  It is important to be aware that there are certain conditions under which this permutation-based FDR is not stable. This includes the case where all observed tests achieve the significance level and the case where there are zero positive tests among the permuted results. In the latter case, we take the conservative approach of setting the number of positive permutation tests to one. Both these cases produce an NA for the upper confidence limit. The latter problem may be alleviated by increasing the number of permutations conducted.
}
\value{
  A dataframe which includes the following columns:
  \item{p.cit }{CIT (omnibus intersection/union test, IUT) p-value, max of component p-values}
  \item{q.cit }{CIT (omnibus FDR intersection/union test, IUT) q-value: permutation-based q-value for p.cit, combining the component FDRs evaluated at the common threshold p.cit (see Details)}
  \item{q.cit.ll }{Lower 95 percent confidence limit (if cl=.95) for q.cit estimate}
  \item{q.cit.ul }{Upper 95 percent confidence limit (if cl=.95) for q.cit estimate}
  \item{q.TaL }{component q-value (FDR) for the test of association between T and L}
  \item{q.ll.TaL }{Lower 95 percent confidence limit (if cl=.95) for q.TaL estimate}
  \item{q.ul.TaL }{Upper 95 percent confidence limit (if cl=.95) for q.TaL estimate}
  \item{q.TaGgvL }{component q-value (FDR) for the test of association between T and G|L}
  \item{q.ll.TaGgvL }{Lower 95 percent confidence limit (if cl=.95) for q.TaGgvL estimate}
  \item{q.ul.TaGgvL }{Upper 95 percent confidence limit (if cl=.95) for q.TaGgvL estimate}
  \item{q.GaLgvT }{component q-value (FDR) for the test of association between G and L|T (for cit.cox results, G and L given the event indicator; absent with test3 = FALSE)}
  \item{q.ll.GaLgvT }{Lower 95 percent confidence limit (if cl=.95) for q.GaLgvT estimate}
  \item{q.ul.GaLgvT }{Upper 95 percent confidence limit (if cl=.95) for q.GaLgvT estimate}
  \item{q.LiTgvG }{component q-value (FDR) for test 4, that the association between L and T given G is weaker than expected without mediation}
  \item{q.ll.LiTgvG }{Lower 95 percent confidence limit (if cl=.95) for q.LiTgvG estimate}
  \item{q.ul.LiTgvG }{Upper 95 percent confidence limit (if cl=.95) for q.LiTgvG estimate}
  The component p-values are also returned, in columns named as in the input.
}
\references{
 Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. btw135. PMID: 27153715.
 Millstein J, Volfson D. 2013. Computationally efficient permutation-based confidence interval estimation for tail-area FDR. Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.

 Efron B, Tibshirani R, Storey JD, Tusher V. 2001. Empirical Bayes analysis of a microarray experiment. Journal of the American Statistical Association 96(456):1151-1160.

 Efron B. 2008. Microarrays, empirical Bayes and the two-groups model. Statistical Science 23(1):1-22.

 Storey JD. 2002. A direct approach to false discovery rates. Journal of the Royal Statistical Society, Series B 64(3):479-498.

 Storey JD. 2003. The positive false discovery rate: a Bayesian interpretation and the q-value. Annals of Statistics 31(6):2013-2035.
}
\author{
  Joshua Millstein
}

\examples{
# Sample Size
ss = 100
n.perm = 20
perm.index = matrix(NA, nrow=ss, ncol=n.perm )
for( j in 1:ncol(perm.index) ) perm.index[, j] = sample( 1:ss )

n.tests = 20
myresults = vector('list', n.tests)

for( tst in 1:n.tests ){

	# Errors
	e1 = matrix(rnorm(ss),ncol=1)
	e2 = matrix(rnorm(ss),ncol=1)

	# Simulate genotypes, gene expression, and clinical traits
	L = matrix(rbinom(ss,2,.5),ncol=1)
	G =  matrix(.5*L + e1,ncol=1)
	T =  matrix(.3*G + e2,ncol=1)
	T = ifelse( T > median(T), 1, 0 )

	myresults[[ tst ]] = cit.bp(L, G, T, perm.index=perm.index, n.perm=n.perm)
}
fdr.cit( myresults )

}

\keyword{ nonparametric }
