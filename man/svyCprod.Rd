\name{svyCprod}
\alias{svyCprod}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Computations for survey variances}
\description{
Computes the sum of products needed for the variance of survey sample estimators.    
}
\usage{
svyCprod(x, strata, psu)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A vector or matrix}
  \item{strata}{A vector of stratum indicators, or \code{NULL}}
  \item{psu}{A vector of cluster indicators or \code{NULL}}
}
\details{
The data are centered within each stratum.  The observations for each cluster are added, and the outer product is taken of the row vector resulting for each cluster.  This is added within strata, multiplied by a degrees-of-freedom correction, and added across strata.
}
\value{
  A matrix
}
\author{Thomas Lumley}


\seealso{\code{\link{svydesign}}, \code{\link{svy.varcoef}} }

\examples{
}
\keyword{utilities}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
