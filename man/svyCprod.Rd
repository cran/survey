\name{svyCprod}
\alias{svyCprod}
\alias{onestage}
\alias{onestrat}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Computations for survey variances}
\description{
Computes the sum of products needed for the variance of survey sample
estimators.  \code{svyCprod} is used for survey design objects from
before version 2.9, \code{onestage} is called by \code{\link{svyrecvar}}
for post-2.9 design objects.
}
\usage{
svyCprod(x, strata, psu, fpc, nPSU,certainty=NULL, postStrata=NULL,
      lonely.psu=getOption("survey.lonely.psu"))
onestage(x, strata, clusters, nPSU, fpc,
      lonely.psu=getOption("survey.lonely.psu"),stage=0,cal)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A vector or matrix}
  \item{strata}{A vector of stratum indicators (may be \code{NULL} for \code{svyCprod})}
  \item{psu}{A vector of cluster indicators (may be \code{NULL})}
  \item{clusters}{A vector of cluster indicators }
  \item{fpc}{A data frame (\code{svyCprod}) or vector (\code{onestage})
    of population stratum sizes, or \code{NULL}}
  \item{nPSU}{Table (\code{svyprod}) or vector (\code{onestage})
    of original sample stratum sizes (or \code{NULL})}
  \item{certainty}{logical vector with stratum names as names. If
    \code{TRUE} and that stratum has a single PSU it is a certainty PSU}
  \item{postStrata}{Post-stratification variables}
  \item{lonely.psu}{One of \code{"remove"}, \code{"adjust"},
    \code{"fail"}, \code{"certainty"}, \code{"average"}. See Details
    below}
  \item{stage}{Used internally to track the depth of recursion}
  \item{cal}{Used to pass calibration information at stages below the population}
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
 
 See \code{\link{postStratify}} for information about
 post-stratification adjustments.
 
The variance formula gives 0/0 if a stratum contains only one sampling
unit. If the \code{certainty} argument specifies that this is a PSU
sampled with probability 1 (a "certainty" PSU) then it does not
contribute to the variance (this is correct only when there is no 
subsampling within the PSU -- otherwise it should be defined as a 
pseudo-stratum).  If \code{certainty} is \code{FALSE} for
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
\code{svydesign} will detect certainty PSUs, making this option
unnecessary.

For strata with a single PSU in a subset (domain) the variance formula
gives a value that is well-defined and positive, but not typically
correct. If \code{options("survey.adjust.domain.lonely")} is \code{TRUE}
and \code{options("survey.lonely.psu")} is \code{"adjust"} or
\code{"average"}, and no post-stratification or G-calibration has been
done, strata with a single PSU in a subset will be treated like those
with a single PSU in the sample.  I am not aware of any theoretical
study of this procedure, but it should at least be conservative.

}
\value{
  A covariance matrix
}
\author{Thomas Lumley}

\references{Binder, David A. (1983).  On the variances of asymptotically normal estimators from complex surveys.  International Statistical Review, 51, 279- 292. }

\seealso{\code{\link{svydesign}}, \code{\link{svyrecvar}}, \code{\link{surveyoptions}}, \code{\link{postStratify}} }

\keyword{utilities}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
