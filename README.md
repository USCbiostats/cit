# cit: Causal Inference Test

The causal inference test (CIT) is a hypothesis test for mediation: whether a
potential mediator G (for example gene expression or DNA methylation) carries
part of a known causal effect of an instrumental variable L (for example a DNA
variant) on an outcome T. It returns a p-value, the maximum of the p-values for
its component conditions (an intersection-union test), and permutation-based FDR
values with confidence intervals when many tests are run.

* `cit.cp()` continuous outcome, `cit.bp()` binary outcome, `cit.cox()` survival outcome
* one mediator or several tested jointly; covariates for the outcome (`CT`) and mediator (`CG`) models
* `robust = TRUE` (CIT v2) prevents a pair from being called causal in both directions;
  recommended when reverse causation is a concern
* `fdr.cit()` for permutation-based q-values across many tests

## Installation

```r
# CRAN release
install.packages("cit")
# development version
remotes::install_github("USCbiostats/cit")
```

## References

Millstein J, Chen GK, Breton CV (2016). cit: hypothesis testing software for mediation
analysis in genomic applications. *Bioinformatics* 32(15):2364-2365.

Millstein J, Zhang B, Zhu J, Schadt EE (2009). Disentangling molecular relationships with
a causal inference test. *BMC Genetics* 10:23.
