\name{cit-package}
\alias{cit-package}
\alias{cit}
\docType{package}
\title{
Causal Inference Test
}
\description{
This package implements a formal statistical hypothesis test for causal mediation. For example, it could be used to test for mediation of a known causal association between a DNA variant, the 'instrumental variable', and a clinical outcome or phenotype by gene expression or DNA methylation, the potential mediator. Another example would be testing mediation of the effect of a drug on a clinical outcome by the molecular target. The hypothesis test generates a p-value or permutation-based false discovery rate (FDR) value with confidence intervals to quantify uncertainty in the causal inference. The outcome can be continuous (\code{\link{cit.cp}}), binary (\code{\link{cit.bp}}) or a censored survival time (\code{\link{cit.cox}}); the potential mediator is continuous and may be multivariate; and the instrumental variable can be continuous or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
}
\details{
This package implements a statistical framework in which existing notions of causal mediation are formalized into a hypothesis test, thus providing a standard quantitative measure of uncertainty in the form of a p-value and optionally, permutation-based FDR. We treat the causal inference as a 'chain' of mathematical conditions that must be satisfied to conclude that the potential mediator is causal for the trait, where the inference is only as good as the weakest link in the chain. P-values are computed for the component conditions. The intersection-union test, in which a series of statistical tests are combined to form an omnibus test, is then employed to generate the overall test result as the maximum of the component p-values. If we let L denote a locus of interest or other instrumental variable, T denote an outcome trait, and G denote a potential mediator, then the four component conditions are 1) L and T are associated, 2) G and T are associated given L, 3) L and G are associated given T, and 4) the association between L and T given G is weaker than expected if G were not a mediator. Condition 4 is tested against null data generated without mediation (Millstein et al. 2009), so the test addresses mediation and can detect it when L also has a direct effect on T.

Version 3.0.0 adds survival outcomes, multiple mediators, separate covariates for the outcome and mediator models, CIT v2 (\code{robust = TRUE}), which prevents a pair of variables from being called causal in both directions and is recommended when reverse causation is a concern, and a choice between empirical and parametric versions of test 4.
}
\author{
Joshua Millstein, Mingzhi Ye

Maintainer: Joshua Millstein <joshua.millstein@usc.edu>
}
\references{
Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. btw135. PMID: 27153715.

Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
}
\keyword{ htest }
\keyword{ nonparametric }
