\name{cit.cox}
\alias{cit.cox}
\title{
  Causal Inference Test for a Survival Outcome
}
\description{
   This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a censored survival time. If the number of permutations is greater than zero, then the results can be used with \code{\link{fdr.cit}} to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is a possibly censored survival time, the potential mediator is continuous and may be multivariate, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
}
\usage{
cit.cox(L, G, T, E, CT = NULL, n.resampl = 50, n.perm = 0, perm.index = NULL,
        rseed = NULL, CG = NULL, test4 = c("empirical", "parametric"),
        maxit = 10000, test3 = TRUE)
}
\arguments{
  \item{L}{Vector or n x p design matrix representing the instrumental variable(s).}
  \item{G}{Continuous vector, or n x k matrix of k continuous potential mediators tested jointly.}
  \item{T}{Numeric vector of follow-up times.}
  \item{E}{Event indicator, coded 1 for an event and 0 for censoring.}
  \item{CT}{Vector or n x p matrix of adjustment covariates for the models with \code{T} as the outcome (tests 1, 2 and 4).}
  \item{n.resampl}{Number of null instances of the statistic for the last test used to estimate the null distribution when \code{test4 = "parametric"} (Millstein et al. 2009). Ignored for the empirical test.}
  \item{n.perm}{If greater than 0, \code{n.perm} permutations of each component test are conducted, generating results under the null for use with \code{\link{fdr.cit}}.}
  \item{perm.index}{Optional n x n.perm matrix whose columns are permutations of \code{1:n}. To account for dependencies among tests when computing FDR confidence intervals, the same permutations must be applied to every CIT in the analysis; pass the same \code{perm.index} to each call. Rows refer to all n observations, including any with missing values. When supplied, the number of permutations is \code{ncol(perm.index)} and \code{n.perm} is not used.}
  \item{rseed}{Optional integer seed. When supplied, results are reproducible and the caller's random number stream is left unchanged. When \code{NULL}, the current stream is used, so \code{set.seed()} before the call also makes results reproducible.}
  \item{CG}{Vector or n x p matrix of adjustment covariates for the models with \code{G} as the outcome (the mediator model used to generate null data for the last test).}
  \item{test4}{How the p-value for the last test is obtained from null data: \code{"empirical"} (default) counts null statistics; \code{"parametric"} fits a non-central chi-square distribution to the likelihood-ratio statistics (Millstein et al. 2009).}
  \item{maxit}{Maximum number of null statistics for the empirical version of the last test.}
  \item{test3}{If \code{TRUE} (default), the CIT includes the test of G given the event indicator, which guards against reverse causation (see Details). \code{FALSE} omits it.}
}
\details{
Let L denote the instrumental variable, G the potential mediator and T the outcome. For a survival outcome the CIT is an intersection-union test of four component conditions, each with its own p-value:
\enumerate{
  \item T is associated with L (\code{p_TassocL}; likelihood-ratio test for L in the Cox proportional hazards model T ~ CT + L, Efron ties).
  \item T is associated with G given L (\code{p_TassocGgvnL}; likelihood-ratio test for G in the Cox model T ~ CT + L + G).
  \item G is associated with L given the event indicator E (\code{p_GassocLgvnT}; F test for L in the linear model G ~ E + CG + L).
  \item The association between L and T given G is weaker than expected if G were not a mediator (\code{p_LindTgvnG}).
}
The third test guards against reverse causation, here the event changing G, as when tissue is collected at or after the event. If G depends on L only through the event, L and G are independent given E and the test does not reject. If G is measured before any event can occur, reverse causation is ruled out by design, and the test can be omitted with \code{test3 = FALSE}; \code{p_cit} is then the maximum of the other three components.

The omnibus p-value, \code{p_cit}, is the maximum of the component p-values, representing the probability of the data if at least one of the component null hypotheses is true.

The last test compares the observed likelihood-ratio statistic for L in the Cox model T ~ CT + G + L with statistics computed from null data in which G is replaced by G* = fitted(G ~ L + CG) + permuted residuals. G* keeps its association with L but has no link to T except through L, so these null data represent no mediation. The test is therefore a test for mediation, not for the absence of a direct effect of L on T: it can reject when part of the effect of L on T is mediated by G and part is direct. With \code{test4 = "empirical"}, the p-value is (r + 1) / (B + 1), where r of B null statistics are at or below the observed statistic; resampling begins with 1000 null statistics and continues, up to \code{maxit}, only while this p-value could determine \code{p_cit}. With \code{test4 = "parametric"}, \code{n.resampl} null statistics give a method-of-moments estimate of the non-centrality, lambda = mean(X*) - df, for a non-central chi-square distribution; each statistic is mapped to a normal score through the fitted distribution, and the p-value is the normal probability of the observed score given the mean and standard deviation of the null scores.

Each component test uses the observations complete for its own variables: the first test uses L, T, E and CT; the second L, G, T, E and CT; the third L, G, T, E and CG; the last all of L, G, T, E, CT and CG.

If permutations are conducted, the results are returned as a data frame in which each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed data. Each permutation is a single reordering of the n observations, applied to L in the first and third tests, to G in the second, and to the mediator-model residuals in the last, so that dependence among the component tests, and among CITs that share a \code{perm.index}, is preserved. In the first three tests the variable is permuted before observations with missing values are removed, so missing values move with their rows; in the last test the permutation is restricted to its observations, keeping each pairing whose source observation is complete. These results can be aggregated across tests and passed to \code{\link{fdr.cit}}.
}
\value{
A data frame. If \code{n.perm = 0} it has one row; otherwise it has a first column \code{perm} (0 for the observed data, then 1, ..., n.perm). Columns:
  \item{p_cit}{CIT (omnibus) p-value.}
  \item{p_TassocL}{component p-value for the test of association between T and L.}
  \item{p_TassocGgvnL}{component p-value for the test of association between T and G given L.}
  \item{p_GassocLgvnT}{component p-value for the test of association between G and L given the event indicator (absent when \code{test3 = FALSE}).}
  \item{p_LindTgvnG}{component p-value for the last test, that the association between L and T given G is weaker than expected without mediation.}
}
\references{
 Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. btw135. PMID: 27153715.

 Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
}
\author{
  Joshua Millstein, Mingzhi Ye
}
\seealso{\code{\link{cit.cp}}, \code{\link{cit.bp}}, \code{\link{fdr.cit}}}
\examples{
set.seed(1)
ss <- 200
L <- matrix(rbinom(ss * 2, 2, .5), ncol = 2)
G <- drop(L \%*\% c(.3, .3)) + rnorm(ss)
T <- rexp(ss, exp(.5 * G))
E <- rbinom(ss, 1, .8)

cit.cox(L, G, T, E, maxit = 2000)
cit.cox(L, G, T, E, test4 = "parametric")
cit.cox(L, G, T, E, test4 = "parametric", test3 = FALSE)  # G measured before any event
}
\keyword{ htest }
\keyword{ survival }
