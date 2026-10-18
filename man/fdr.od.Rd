\name{fdr.od}
\alias{fdr.od}
\title{
  Permutation-Based FDR and Confidence Interval
}
\description{
  This function can be used to estimate FDR, corresponding confidence interval, and pi0, the proportion of true null hypotheses, given a selected significance threshold, and results from permuted data.
}
\usage{
fdr.od(obsp, permp, pnm, ntests, thres, cl=.95, od=NA)
}

\arguments{
  \item{obsp}{
     observed vector of p-values.
}
  \item{permp}{
     list of dataframes that include a column of permutation p-values (or statistics) in each. The length of the list permp = number of permutations.
}
  \item{pnm}{
     name of column in each list component dataframe that includes p-values (or statistics).
}
  \item{ntests}{
     total number of observed tests, which is usually the same as the length of obsp and the number of rows in each permp dataframe. However, this may not be the case if results were filtered by a p-value threshold or statistic threshold. If filtering was conducted then thres must be smaller (more extreme) than the filtering criterion.
}
  \item{thres}{
     significance threshold.
}
  \item{cl}{
     confidence level for FDR confidence interval (default is .95).
}
  \item{od}{
     Over-dispersion parameter. Setting the over-dispersion parameter to one corresponds to the assumption that all tests (referring to multiple omnibus tests conducted using different sets of variables) are independent. The default (c1=NA) is to use an empirically estimated over-dispersion parameter based on the permutation results.
}
}
\details{
  If a very large number of tests are conducted, it may be useful to filter results, that is, save only results of those tests that meet some relaxed nominal significance threshold. This alleviates the need to record results for tests that are clearly non-significant. Results from fdr_od() are valid as long as thres < the relaxed nomimal significance threshold for both observed and permuted results. It is not necessary for the input to fdr_od() to be p-values, however, fdr_od() is designed for statistics in which smaller values are more extreme than larger values as is the case for p-values. Therefore, if raw statistics are used, then a transformation may be necessary to insure that smaller values are more likely associated with false null hypotheses than larger values.
  
  If there are zero positive tests among any of the permuted results (s.perm = 0) at the specified significance threshold (thres), then a conservative approximation is conducted by setting s.perm = 1 for the FDR calculation. This approximation is necessary for the permutation-based FDR to be estimable. This situation can be identified by noting that s.perm = 0 in the output. Increasing the number of permutations conducted may result in positive permutation tests, and thereby overcome the problem. In certain situations, for instance when a large proportion of tests meet the significance threshold, pi0 is estimated to be very small, and thus has a large influence on the FDR estimate. To limit this influence, pi0 is constrained to be .5 or greater, resulting in a more conservative estimate under these conditions. 
}
\value{
  A vector which includes:
  \item{FDR }{FDR point estimate}
  \item{ll }{lower confidence limit}
  \item{ul }{upper confidence limit}
  \item{pi0 }{proportion of true null hypotheses}
  \item{od }{overdispersion parameter}
  \item{s.obs }{observed number of positive tests}
  \item{s.perm }{total number of positive tests summed across all permuted result sets}
}
\references{
  Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. btw135. PMID: 27153715.
  Millstein J, Volfson D. 2013. Computationally efficient permutation-based confidence interval estimation for tail-area FDR. Frontiers in Genetics | Statistical Genetics and Methodology 4(179):1-11.
}
\author{
  Joshua Millstein
}

\examples{
nrow_=100
ncol_=100
X = as.data.frame(matrix(rnorm(nrow_*ncol_),nrow=nrow_,ncol=ncol_))
Y = as.data.frame(matrix(rnorm(nrow_*ncol_),nrow=nrow_,ncol=ncol_))
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
for(p_ in 1:nperm){
	X1 = X[order(runif(ncol_)),]
	perml[[p_]] = myanalysis(X1,Y)
}

## FDR results
fdr.od(obs$pvalue,perml,"pvalue",ncol_,.05)

}

\keyword{ htest }
\keyword{ nonparametric }
