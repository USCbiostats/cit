\name{cit.effects}
\alias{cit.effects}
\title{
  Direct and Indirect Effect Estimates
}
\description{
  Estimates the total, direct and indirect effects of the instrumental variable L on the outcome T through the mediator G, with delta-method confidence intervals, using the regression estimators of VanderWeele and Vansteelandt in the form given by Valeri and VanderWeele (2013). It is intended to accompany \code{\link{cit.cp}}, \code{\link{cit.bp}} and \code{\link{cit.cox}}: the CIT tests whether the data support mediation, and these estimates describe its size.
}
\usage{
cit.effects(L, G, T, E = NULL, CT = NULL, CG = NULL,
            outcome = c("continuous", "binary", "survival"), conf.level = 0.95)
}
\arguments{
  \item{L}{Vector or n x p design matrix representing the instrumental variable(s). Effects are reported for each column.}
  \item{G}{Continuous vector: a single potential mediator.}
  \item{T}{Outcome: a continuous vector, a binary vector coded 0/1, or a vector of follow-up times.}
  \item{E}{Event indicator (1 = event, 0 = censored), for \code{outcome = "survival"} only.}
  \item{CT}{Vector or n x p matrix of adjustment covariates for the outcome model.}
  \item{CG}{Vector or n x p matrix of adjustment covariates for the mediator model.}
  \item{outcome}{Type of outcome: \code{"continuous"} (linear regression, as in \code{cit.cp}), \code{"binary"} (logistic regression, as in \code{cit.bp}) or \code{"survival"} (Cox regression with Efron ties, as in \code{cit.cox}).}
  \item{conf.level}{Confidence level of the intervals.}
}
\details{
Two models are fitted to the observations complete for every variable used: the mediator model G ~ L + CG (linear) and the outcome model T ~ L + G + CT (linear, logistic or Cox). As in the CIT, the models have no L x G interaction. Let b1 be the coefficient of a column of L in the mediator model, and t1 and t2 the coefficients of that column and of G in the outcome model. For a one-unit increase in that column, with the other columns held fixed:
\itemize{
  \item direct effect = t1;
  \item indirect effect = t2 b1;
  \item total effect = t1 + t2 b1.
}
For a continuous outcome these are mean differences. For a binary outcome they are log odds ratios and for a survival outcome log hazard ratios, reported as odds ratios or hazard ratios. The ratio-scale decomposition assumes a rare outcome (VanderWeele and Vansteelandt 2010; VanderWeele 2011), and a warning is given when more than 10\% of observations have the event. The proportion mediated is indirect / total for a continuous outcome and OR_D (OR_I - 1) / (OR_D OR_I - 1) on the ratio scales (with HR for OR).

Standard errors are by the delta method, treating the two models as independent (Valeri and VanderWeele 2013). The p-values are Wald tests of a zero effect (on the log scale for ratios). The p-value for the indirect effect is a Sobel-type test; it is not the CIT, which also tests the direction of the relationship.

These estimates rely on stronger assumptions than the CIT. Interpreting them as natural direct and indirect effects requires no unmeasured confounding of the L-T, L-G and G-T relationships given the covariates, and no confounder of G and T that is itself affected by L. When L is a genetic variant, its relationships with G and T are protected by random assortment, but the G-T relationship is not. Measurement error in G biases the indirect effect towards zero and the direct effect away from zero.
}
\value{
A data frame with one row per column of L and effect: \code{L} (column name), \code{effect} (\code{"total"}, \code{"direct"}, \code{"indirect"}, \code{"prop.mediated"}), \code{scale} (\code{"difference"}, \code{"OR"}, \code{"HR"} or \code{"proportion"}), \code{estimate}, \code{lower}, \code{upper}, \code{se} (on the log scale for ratios), \code{p} and \code{n} (observations used).
}
\references{
  VanderWeele TJ, Vansteelandt S. 2009. Conceptual issues concerning mediation, interventions and composition. Statistics and Its Interface 2:457-468.

  VanderWeele TJ, Vansteelandt S. 2010. Odds ratios for mediation analysis for a dichotomous outcome. American Journal of Epidemiology 172:1339-1348.

  VanderWeele TJ. 2011. Causal mediation analysis with survival data. Epidemiology 22:582-585.

  Valeri L, VanderWeele TJ. 2013. Mediation analysis allowing for exposure-mediator interactions and causal interpretation: theoretical assumptions and implementation with SAS and SPSS macros. Psychological Methods 18:137-150.
}
\seealso{
  \code{\link{cit.cp}}, \code{\link{cit.bp}}, \code{\link{cit.cox}}
}
\examples{
set.seed(1)
n <- 500
L <- rbinom(n, 2, 0.3)
G <- 0.5 * L + rnorm(n)
T <- 0.3 * G + 0.1 * L + rnorm(n)
cit.effects(L, G, T)
}
\keyword{ htest }
\keyword{ regression }
