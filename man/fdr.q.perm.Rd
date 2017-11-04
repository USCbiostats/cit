\name{fdr.q.perm}
\alias{fdr.q.perm}
\title{
  Nonparametric permutation-based tail-area FDR Values, q-values
}
\description{
   Given a vector of p-values, this function uses fdr.od() to compute tail-area FDR to assign q-values, minimum tail-area FDR, for each p-value. This function is used by the function fdr.cit to generate q-values for the three component tests, TassocL, TassocGgvnL, and GassocLgvnT.
}
\usage{
   fdr.q.perm( obs.p, perml, pname, ntests, cl=.95, od=NA )
}

\arguments{
  \item{obs.p}{
     observed vector of p-values.
}
  \item{perml}{
     list of dataframes that include a column of permutation p-values (or statistics) in each. The length of the list permp = number of permutations.
}
  \item{pname}{
     name of column in each list component dataframe that includes p-values (or statistics).
}
  \item{ntests}{
     total number of observed tests, which is usually the same as the length of obsp and the number of rows in each permp dataframe. However, this may not be the case if results were filtered by a p-value threshold or statistic threshold. If filtering was conducted then thres must be smaller (more extreme) than the filtering criterion.
}
  \item{cl}{
     confidence level for FDR confidence interval (default is .95).
}
  \item{od}{
     Over-dispersion parameter. Setting the over-dispersion parameter to one corresponds to the assumption that all tests (referring to multiple omnibus tests conducted using different sets of variables) are independent. The default (c1=NA) is to use an empirically estimated over-dispersion parameter based on the permutation results.
}
}
\details{
  Following Storey and Tibshirani (2003), the q-value is set to the minimum FDR for p-values greater than or equal to the observed p-value.
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
rowno=100
colno=100
X = as.data.frame(matrix(rnorm(rowno*colno),nrow=rowno,ncol=colno))
Y = as.data.frame(matrix(rnorm(rowno*colno),nrow=rowno,ncol=colno))
nperm = 10

myanalysis = function(X,Y){
	ntests = ncol(X)
	rslts = as.data.frame(matrix(NA,nrow=ntests,ncol=2))
	names(rslts) = c("ID","pvalue")
	rslts[,"ID"] = 1:ntests
	for(i in 1:ntests){
		fit = cor.test(X[,i],Y[,i],na.action="na.exclude",
			alternative="two.sided",method="pearson")
		rslts[i,"pvalue"] = fit$p.value
	}
	return(rslts)
} # End myanalysis

# Generate observed results
obs = myanalysis(X,Y)

## Generate permuted results
perml = vector('list',nperm)
for(perm in 1:nperm){
	X1 = X[order(runif(colno)),]
	perml[[perm]] = myanalysis(X1,Y)
}

## FDR results
fdr.q.perm(obs$pvalue,perml,"pvalue",colno)

}

\keyword{ q-value }
