\name{bootweights}
\alias{bootweights}
\alias{bootstratum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute survey bootstrap weights }
\description{
Bootstrap weights for infinite populations are created by sampling with
replacement from the PSUs in each stratum.  For finite populations an
empirical population is created as described by Canty and Davision
(1999).  The bootstrap is primarily useful when most strata have many
PSUs.

This bootstrap is strictly appropriate only when the first stage of
sampling is a simple or stratified random sample of PSUs with or without
replacement, and not (eg) for PPS sampling.
}
\usage{
bootweights(strata, psu, replicates = 50, fpc = NULL,
         fpctype = c("population", "fraction", "correction"),
         compress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strata}{Identifier for sampling strata (top level)}
  \item{psu}{Identifier for primary sampling units}
  \item{replicates}{Number of bootstrap replicates}
  \item{fpc}{Finite population correction }
  \item{fpctype}{Is \code{fpc} the population size, sampling fraction,
    or 1-sampling fraction?}
  \item{compress}{Should the replicate weights be compressed?}
}

\value{
 A set of replicate weights
}

\references{Canty AJ, Davison AC. (1999) Resampling-based variance
 estimation for labour force surveys. The Statistician 48:379-391 }

\seealso{\code{\link{as.svrepdesign}}}
\examples{

}
\keyword{survey}% at least one, from doc/KEYWORDS

