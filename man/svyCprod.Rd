\name{svyCprod}
\alias{svyCprod}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Computations for survey variances}
\description{
Computes the sum of products needed for the variance of survey sample estimators.    
}
\usage{
svyCprod(x, strata, psu, fpc, nPSU,certainty=NULL, postStrata=NULL
      lonely.psu=getOption("survey.lonely.psu"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A vector or matrix}
  \item{strata}{A vector of stratum indicators, or \code{NULL}}
  \item{psu}{A vector of cluster indicators or \code{NULL}}
  \item{fpc}{A data frame of population stratum sizes or \code{NULL}}
  \item{nPSU}{Table of original sample stratum sizes (or \code{NULL})}
  \item{certainty}{logical vector with stratum names as names. If
    \code{TRUE} and that stratum has a single PSU it is a certainty PSU}
  \item{postStrata}{Post-stratification variables}
  \item{lonely.psu}{One of \code{"remove"}, \code{"adjust"},
    \code{"fail"}, \code{"certainty"}, \code{"average"}. See Details below}

}
\details{
 The observations for each cluster are added, then centered within each
 stratum and the outer product is taken of the row vector resulting for
 each cluster.  This is added within strata, multiplied by a
 degrees-of-freedom correction and by a finite population correction (if
 supplied) and added across strata.  

 If there are fewer clusters (PSUs) in a stratum than in the original
 design extra rows of zeroes are added to \code{x} to allow the correct
 subpopulation variance to be computed.

The variance formula gives 0/0 if a stratum contains only one sampling
unit. If the \code{certainty} argument specifies that this is a PSU
sampled with probability 1 (a "certainty" PSU) then it does not
contribute to the variance.  If \code{certainty} is \code{FALSE} for
this stratum or is not supplied the result depends on \code{lonely.psu}.

The options are \code{"fail"} to give an error, \code{"remove"} or
\code{"certainty"} to give a variance contribution of 0 for the stratum,
\code{"adjust"} to center the stratum at the grand mean rather than the
stratum mean, and \code{"average"} to assign strata with one PSU the
average variance contribution from strata with more than one PSU.  The
choice is controlled by setting \code{options(survey.lonely.psu)}. If
this is not done the factory default is \code{"fail"}. Using
\code{"adjust"} is conservative, and it would often be better to combine
strata in some intelligent way. The properties of \code{"average"} have
not been investigated thoroughly, but it may be useful when the lonely
PSUs are due to a few strata having PSUs missing completely at random.

The \code{"remove"}and \code{"certainty"} options give the same result,
but \code{"certainty"} is intended for situations where there is only
one PSU in the population stratum, which is sampled with certainty (also
called `self-representing' PSUs or strata). With \code{"certainty"} no
warning is generated for strata with only one PSU.  Ordinarily,
\code{svydesign} will detect certainty PSUs, making this option unnecessary.


When a subset of a survey design has only one PSU in a stratum the
standard formulas give zero for the contribution of that stratum to the
variance.  \code{svyCprod} will warn that this has happened, and if
\code{lonely.psu="adjust"} or \code{lonely.psu="average"} will use the
same adjustment as if the whole survey had only one PSU in that stratum
by design.

}
\value{
  A covariance matrix
}
\author{Thomas Lumley}


\seealso{\code{\link{svydesign}}, \code{\link{svy.varcoef}} }

\examples{
}
\keyword{utilities}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
