\name{cit.bp}
\alias{cit.bp}
\title{
  Causal Inference Test for a Binary Outcome
}
\description{
   This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a trait or clinical outcome. If the number of permutations is greater than zero,  then the results can be used with fdr.cit to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is binary, the potential mediator is continuous, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
}
\usage{
cit.bp( L, G, T, C=NULL, maxit=10000, n.perm=0, perm.index=NULL, rseed=NULL )
}
\arguments{
  \item{L}{vector or nxp design matrix representing the instrumental variable(s). 
}
  \item{G}{continuous vector representing the potential causal mediator.}
  \item{T}{
     Continuous vector representing the clinical trait or outcome of interest.
}
  \item{C}{
     Vector or nxp design matrix representing adjustment covariates. 
}
  \item{maxit}{
      Maximum number of iterations to be conducted for the conditional independence test, test 4, which is permutation-based. The minimum number of permutations conducted is 1000, regardless of maxit. Increasing maxit will increase the precision of the p-value for test 4 if the p-value is small.
}
  \item{n.perm}{
      If n.perm is set to an integer greater than 0, then n.perm permutations for each component test will be conducted (randomly permuting the data to generate results under the null). 
}
  \item{perm.index}{
      This item is only important when the CIT is conducted multiple times, leading to a multiple testing issue, that is addressed by computing FDR using the fdr.cit() function. To accurately account for dependencies among tests when computing FDR confidence intervals, the permutations must be the same for all tests. Thus, if 100 permutations are conducted for 500 CIT test scenarios, each of these 100 permutations are applied to all 500 tests. This is achieved by passing an argument, perm.index, which is an n row by n.perm column dataframe or matrix of permutation indices, where n is the number of observations and n.perm the number of permutations. Each column of perm.index includes a random permutation of 1:n, and this same perm.index object would passed to all 500 CIT tests. fdr.cit() would then be used to compute q-values  (FDR) and q-value confidence intervals for each test.
}
  \item{rseed}{
      If n.perm > 0, and multiple tests (CITs) are being conducted, setting rseed to the same integer for all tests insures that the permutations will be the same across CITs. This is important for maintaining the observed dependencies among tests for permuted data in order to compute accurate confidence intervals for FDR estimates. 
}
}
\details{
  The omnibus p-value, p_cit, is the maximum of the component p-values, an intersection-union test, representing the probability of the data if at least one of the component null hypotheses is true. For component test 4, rather than using the semiparametric approach proposed by Millstein et al. (2009), here it is estimated completely by permutation, resulting in an exact test. If permutations are conducted by setting n.perm to a value greater than zero, then the results are provided in matrix (dataframe) form where each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed or non-permuted analysis. These results can then be aggregated across multiple cit.bp tests and input to the function fdr.cit to generate component test FDR values (q-values) as well as omnibus q-values with confidence intervals that correspond to the p_cit omnibus p-values.
}
\value{
  A dataframe which includes the following columns:
  \item{ perm }{Indicator for permutation results. Zero indicates that the data were not permuted and subsequent rows include an integer greater than zero for each permutation conducted.}
  \item{p_cit }{CIT (omnibus) p-value}
  \item{p_TassocL }{component p-value for the test of association between T and L.}
  \item{p_TassocGgvnL }{component p-value for the test of association between T and G|L.}
  \item{p_GassocLgvnT }{component p-value for the test of association between G and L|T.}
  \item{p_LindTgvnG }{component p-value for the equivalence test of L ind T|G}
}
\references{
 Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. btw135. PMID: 27153715.
 Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
}
\author{
  Joshua Millstein
}

\examples{
# Sample Size
ss = 100

# Errors
e1 = matrix(rnorm(ss),ncol=1)
e2 = matrix(rnorm(ss),ncol=1)

# Simulate genotypes, gene expression, covariates and a clinical trait
L = matrix(rbinom(ss*3,2,.5),ncol=3)
G =  matrix( apply(.4*L, 1, sum) + e1,ncol=1)
T =  matrix(.3*G + e2,ncol=1)
T = ifelse( T > median(T), 1, 0 )
C =  matrix(matrix(rnorm(ss*2),ncol=1),ncol=2)

n.perm = 5
perm.index = matrix(NA, nrow=ss, ncol=n.perm )
for( j in 1:ncol(perm.index) ) perm.index[, j] = sample( 1:ss )

results = cit.bp(L, G, T, perm.index=perm.index, n.perm=n.perm)
results

results = cit.bp(L, G, T)
results

results = cit.bp(L, G, T, C, n.perm=5)
results

results = cit.bp(L, G, T, C)
results
}
\keyword{ nonparametric }
