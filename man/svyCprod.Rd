\name{svyCprod}
\alias{svyCprod}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Computations for survey variances}
\description{
Computes the sum of products needed for the variance of survey sample estimators.    
}
\usage{
svyCprod(x, strata, psu, fpc, lonely.psu=getOption("survey.lonely.psu"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A vector or matrix}
  \item{strata}{A vector of stratum indicators, or \code{NULL}}
  \item{psu}{A vector of cluster indicators or \code{NULL}}
  \item{fpc}{A vector of population stratum sizes or \code{NULL}}
  \item{lonely.psu}{One of \code{"remove"}, \code{"adjust"},
    \code{"fail"}. See Details below}
}
\details{
 The observations for each cluster are added, then centred within each
 stratum and the outer product is taken of the row vector resulting for
 each cluster.  This is added within strata, multiplied by a
 degrees-of-freedom correction and by a finite population correction (if
 supplied) and added across strata.

The variance formula gives 0/0 if a stratum contains only one sampling
unit. The options to handle this are \code{"fail"} to give an error,
\code{"remove"} to give a variance contribution of 0 for the stratum,
and \code{"adjust"} to center the stratum at the grand mean rather than
the stratum mean.  The choice is controlled by setting
\code{options(survey.lonely.psu)}. If this is not done the factory
default is \code{"fail"}. Using \code{"adjust"} is conservative, and it
would often be better to combine strata in some intelligent way. The
\code{"remove"} option is included because it might be needed for
compatibility.  }
\value{
  A covariance matrix
}
\author{Thomas Lumley}


\seealso{\code{\link{svydesign}}, \code{\link{svy.varcoef}} }

\examples{
}
\keyword{utilities}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
