\name{cit.bp}
\alias{cit.bp}
\title{
  Causal Inference Test for a Binary Outcome
}
\description{
   This function implements a formal statistical hypothesis test, resulting in a p-value, to quantify uncertainty in a causal inference pertaining to a measured factor, e.g. a molecular species, which potentially mediates a known causal association between a locus or other instrumental variable and a binary trait. If the number of permutations is greater than zero, then the results can be used with \code{\link{fdr.cit}} to generate permutation-based FDR values (q-values) that are returned with confidence intervals to quantify uncertainty in the estimate. The outcome is binary, the potential mediator is continuous and may be multivariate, and the instrumental variable can be continuous, discrete (such as coding a SNP 0, 1, 2), or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
}
\usage{
cit.bp(L, G, T, CT = NULL, maxit = 10000, n.perm = 0, perm.index = NULL,
       rseed = NULL, CG = NULL, robust = FALSE,
       test4 = c("empirical", "parametric"), n.resampl = 50, C = NULL)
}
\arguments{
  \item{L}{Vector or n x p design matrix representing the instrumental variable(s).}
  \item{G}{Continuous vector, or n x k matrix of k continuous potential mediators tested jointly.}
  \item{T}{Binary vector, coded 0/1, representing the clinical trait or outcome of interest.}
  \item{CT}{Vector or n x p matrix of adjustment covariates for the models with \code{T} as the outcome (tests 1, 2 and 4). In cit 2.3.2 and earlier this argument was called \code{C}.}
  \item{n.resampl}{Number of null instances of the test 4 statistic used to estimate the null distribution when \code{test4 = "parametric"} (Millstein et al. 2009). Ignored for the empirical test.}
  \item{n.perm}{If greater than 0, \code{n.perm} permutations of each component test are conducted, generating results under the null for use with \code{\link{fdr.cit}}.}
  \item{perm.index}{Optional n x n.perm matrix whose columns are permutations of \code{1:n}. To account for dependencies among tests when computing FDR confidence intervals, the same permutations must be applied to every CIT in the analysis; pass the same \code{perm.index} to each call. Rows refer to all n observations, including any with missing values.}
  \item{rseed}{Optional integer seed. When supplied, results are reproducible and the caller's random number stream is left unchanged. When \code{NULL}, the current stream is used, so \code{set.seed()} before the call also makes results reproducible.}
  \item{CG}{Vector or n x p matrix of adjustment covariates for the models with \code{G} as the outcome (test 3 and the mediator model used to generate test 4 null data).}
  \item{robust}{If \code{TRUE}, use CIT v2, in which test 3 is referred to a non-central F distribution (see Details). Recommended when reverse causation (T affecting G) is a concern. Available for a single mediator.}
  \item{test4}{How the test 4 p-value is obtained from null data: \code{"empirical"} (default) counts null statistics; \code{"parametric"} fits a non-central F distribution (Millstein et al. 2009).}
  \item{maxit}{Maximum number of null statistics for the empirical test 4.}
  \item{C}{Deprecated name for \code{CT}.}
}
\details{
Let L denote the instrumental variable, G the potential mediator and T the outcome. The CIT is an intersection-union test of four component conditions, each with its own p-value:
\enumerate{
  \item T is associated with L (\code{p_TassocL}; likelihood-ratio test of L in the logistic regression T ~ CT + L).
  \item T is associated with G given L (\code{p_TassocGgvnL}; likelihood-ratio test of G in T ~ CT + L + G).
  \item G is associated with L given T (\code{p_GassocLgvnT}; F test of L in G ~ T + CG + L, or Wilks' lambda when G has several columns).
  \item The association between L and T given G is weaker than expected if G were not a mediator (\code{p_LindTgvnG}).
}
The omnibus p-value, \code{p_cit}, is the maximum of the component p-values, representing the probability of the data if at least one of the component null hypotheses is true.

Test 4 compares the observed likelihood-ratio statistic for L in the logistic regression T ~ CT + G + L with statistics computed from null data in which G is replaced by G* = fitted(G ~ L + CG) + permuted residuals. G* keeps its association with L but has no link to T except through L, so these null data represent no mediation. The test is therefore a test for mediation, not for the absence of a direct effect of L on T: it can reject when part of the effect of L on T is mediated by G and part is direct. With \code{test4 = "empirical"}, the p-value is (r + 1) / (B + 1), where r of B null statistics are at or below the observed statistic; resampling begins with 1000 null statistics and continues, up to \code{maxit}, only while the test 4 p-value could determine \code{p_cit}. With \code{test4 = "parametric"}, \code{n.resampl} null statistics give a method-of-moments estimate of the non-centrality, lambda = mean(X*) - df, for a non-central chi-square distribution; each statistic is mapped to a normal score through the fitted distribution, and the p-value is the normal probability of the observed score given the mean and standard deviation of the null scores.

With \code{robust = TRUE} (CIT v2), the test 3 statistic is referred to a non-central F distribution whose non-centrality is the method-of-moments estimate, lambda = X - df, from the likelihood-ratio statistic X for L in the logistic regression T ~ CT + G + L. For a continuous outcome (\code{\link{cit.cp}}) the corresponding change guarantees that the CIT cannot conclude that G mediates the effect of L on T and also, with the roles of G and T exchanged (and CT and CG exchanged), that T mediates the effect of L on G; for a binary outcome the same construction is used, but the guarantee is approximate.

Each component test uses the observations complete for its own variables: test 1 uses L, T and CT; test 2 L, G, T and CT; test 3 L, G, T and CG; test 4 all of them. With \code{robust = TRUE}, test 3 uses the test 4 observations, so that the test 3 statistic and its non-centrality come from the same data.

If permutations are conducted, the results are returned as a data frame in which each row represents an analysis using a unique permutation, except the first row (perm = 0), which has results from the observed data. Each permutation is a single reordering of the n observations, applied to L in tests 1 and 3 (and, with \code{robust = TRUE}, in the statistic that gives the test 3 non-centrality), to G in test 2, and to the mediator-model residuals in test 4, so that dependence among the component tests, and among CITs that share a \code{perm.index}, is preserved. In tests 1 to 3 the variable is permuted before observations with missing values are removed, so missing values move with their rows; in test 4 the permutation is restricted to the test 4 observations, keeping each pairing whose source observation is complete. These results can be aggregated across tests and passed to \code{\link{fdr.cit}}.
}
\value{
A data frame. If \code{n.perm = 0} it has one row (as in earlier versions); otherwise it has a first column \code{perm} (0 for the observed data, then 1, ..., n.perm). Columns:
  \item{p_cit}{CIT (omnibus) p-value.}
  \item{p_TassocL}{component p-value for the test of association between T and L.}
  \item{p_TassocGgvnL}{component p-value for the test of association between T and G given L.}
  \item{p_GassocLgvnT}{component p-value for the test of association between G and L given T.}
  \item{p_LindTgvnG}{component p-value for test 4, that the association between L and T given G is weaker than expected without mediation.}
}
\references{
 Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. btw135. PMID: 27153715.

 Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
}
\author{
  Joshua Millstein, Mingzhi Ye
}
\seealso{\code{\link{cit.cp}}, \code{\link{cit.cox}}, \code{\link{fdr.cit}}}
\examples{
set.seed(1)
ss <- 200
L <- matrix(rbinom(ss * 3, 2, .5), ncol = 3)
G <- matrix(apply(.3 * L, 1, sum) + rnorm(ss), ncol = 1)
T <- rbinom(ss, 1, plogis(.8 * G))
CT <- matrix(rnorm(ss * 2), ncol = 2)

cit.bp(L, G, T)
cit.bp(L, G, T, CT, robust = TRUE)
cit.bp(L, G, T, CT, n.perm = 5)
}
\keyword{ htest }
\keyword{ nonparametric }
