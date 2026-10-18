\name{iuq}
\alias{iuq}
\title{
  Intersection/Union Q-Value
}
\description{
   Tail-area false discovery rates (FDRs) or 'q-values' are combined across multiple component tests to estimate a single q-value that represents the intersection of alternative hypotheses or the union of null hypotheses. In other words, this is an estimate of the rate at which tests called 'significant' at the specified omnibus q-value include at least one component condition that is called 'significant' due to chance alone.
}
\usage{
iuq( qvec )
}
\arguments{
  \item{qvec}{vector of q-values, each corresponding to a component test. 
}
}
\details{
  The omnibus q-value is greater than or equal to the maximum of the component q-values.
}
\value{
  A single value that represents the rate at which tests called 'significant' at the specified omnibus q-value include at least one component condition that is called 'significant' due to chance alone.
}
\references{
 Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. btw135. PMID: 27153715.
}
\author{
  Joshua Millstein
}

\examples{
# component q-values
qvec = c( .1, .04, .02, .11 )

# omnibus q-value
iuq( qvec )
}
\keyword{ nonparametric }
