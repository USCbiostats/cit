\name{linreg}
\alias{linreg}
\title{
  F Test for Linear Model
}
\description{
  This function is used by cit.cp to compute F test given a continuous outcome and full vs reduced sets of covariates
}
\usage{
linreg( nms.full, nms.redu=NULL, nm.y, mydat )
}
\arguments{
  \item{nms.full}{
     vector of variable names for all covariates included in the full model.
}
  \item{nms.redu}{
     vector of variable names for all covariates included in the reduced model. If nms.redu is NULL then the reduced model is fitted with the intercept only.
}
  \item{nm.y}{
     character string, which is the name of the outcome variable.
}
  \item{mydat}{
     the dataframe that includes all variables with each variable in a column.
}
}
\details{
  An F test is conducted using the glm function by comparing the full and reduced models. This function is called by cit.cp.
}
\value{
  A single p-value is returned.
}
\references{
  Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.
}
\author{
  Joshua Millstein
}

\examples{
ss = 500
cols = 6
nm.y = "y"
nms.full = paste( "x", 1:(cols-1), sep="" )
nms.redu = paste( "x", 1:2, sep="" )

mydat = as.data.frame( matrix( rnorm( ss*cols ), ncol=cols ) )
names( mydat ) = c( nm.y, nms.full )

linreg(nms.full, nms.redu, nm.y, mydat)

}

\keyword{ htest }
\keyword{ nonparametric }
