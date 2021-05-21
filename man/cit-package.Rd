\name{cit-package}
\alias{cit-package}
\docType{package}
\title{
Causal Inference Test
}
\description{
This package implements a formal statistical hypothesis test for causal mediation. For example, it could be used to test for mediation of a known causal association between a DNA variant, the 'instrumental variable', and a clinical outcome or phenotype by gene expression or DNA methylation, the potential mediator. Another example would be testing mediation of the effect of a drug on a clinical outcome by the molecular target. The hypothesis test generates a p-value or permutation-based false discovery rate (FDR) value with confidence intervals to quantify uncertainty in the causal inference. The outcome can be represented by either a continuous or binary variable, the potential mediator is continuous, and the instrumental variable can be continuous or binary and is not limited to a single variable but may be a design matrix representing multiple variables.
}
\details{
\tabular{ll}{
Package: \tab cit\cr
Type: \tab Package\cr
Version: \tab \Sexpr{packageDescription("cit")$Version}\cr
Date: \tab \Sexpr{packageDescription("cit")$Date}\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
}
This package implements a novel statistical framework in which existing notions of causal mediation are formalized into a hypothesis test, thus providing a standard quantitative measure of uncertainty in the form of a p-value and optionally, permutation-based FDR. We treat the causal inference as a 'chain' of mathematical conditions that must be satisfied to conclude that the potential mediator is causal for the trait, where the inference is only as good as the weakest link in the chain. P-values are computed for the component conditions. The Intersection-Union Test, in which a series of statistical tests are combined to form an omnibus test, is then employed to generate the overall test result as the maximum of the component p-values. If we let L denote a locus of interest or other instrumental variable, T denote an outcome trait, and G denote a potential mediator, then the four component conditions are, 1) L and G are associated, 2) L and T are associated, 3) L is associated with G|T, and 4) L is independent of T|G. Test 4 requires an equivalence test and is implemented here using a permutation based approach (as described in Millstein et al. 2009).
}
\author{
Joshua Millstein

Maintainer: Joshua Millstein <joshua.millstein@usc.edu>
Joshua Millstein
}
\references{
Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. btw135. PMID: 27153715.
Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
}

\keyword{ htest }
\keyword{ nonparametric }
