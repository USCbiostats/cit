\name{fdr.q.para}
\alias{fdr.q.para}
\title{
  Parametric tail-area FDR Values, q-values
}
\description{
   Given a vector of p-values, this function will use a parametric expression for tail-area FDR to assign q-values, minimum tail-area FDR, for each p-value. This function is used by the function fdr.cit to generate q-values for the forth test, LindTgvnG.
}
\usage{
   fdr.q.para( pvals )
}

\arguments{
  \item{pvals}{
     Vector of p-values.
}
}
\details{
  The parametric expression is constructed by setting the average number of positive tests among permuted results to p*m, where p is the p-value and m is the number of tests, the length of vector pvals. Following Storey and Tibshirani (2003), the q-value is set to the minimum FDR for p-values greater than or equal to the observed p-value.
}
\value{
  A vector of q-values.
}
\references{
 Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. btw135.
 Millstein J, Volfson D. 2013. Computationally efficient permutation-based confidence interval estimation for tail-area FDR. Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.
 Storey DJ, Tibshirani R. 2003. Statistical significance for genomewide studies. PNAS 100(16):9440-9445.
}
\author{
  Joshua Millstein
}

\examples{
# Sample Size
ss = 100

n.tests = 20
mypvals = rep(NA, n.tests)

for( tst in 1:n.tests ){

	# Errors
	e1 = matrix(rnorm(ss),ncol=1)
	e2 = matrix(rnorm(ss),ncol=1)

	# Simulate genotypes, gene expression, and clinical traits
	L = matrix(rbinom(ss,2,.5),ncol=1)
	G =  matrix(.5*L + e1,ncol=1)
	T =  matrix(.3*G + e2,ncol=1)
	T = ifelse( T > median(T), 1, 0 )

	mypvals[ tst ] = cit.bp(L, G, T)$p_cit
}
fdr.q.para( mypvals )

}

\keyword{ parametric }
