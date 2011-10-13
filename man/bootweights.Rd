\name{bootweights}
\alias{bootweights}
\alias{subbootweights}
\alias{mrbweights}
\alias{bootstratum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute survey bootstrap weights }
\description{
Bootstrap weights for infinite populations ('with replacement' sampling) are created by sampling with
replacement from the PSUs in each stratum. \code{subbootweights()}
samples \code{n-1} PSUs from the \code{n} available (Rao and Wu),
\code{bootweights} samples \code{n} (Canty and Davison). 

For multistage designs or those with large sampling fractions,
\code{mrbweights} implements Preston's multistage rescaled
bootstrap. The multistage rescaled bootstrap is still useful for
single-stage designs with small sampling fractions, where it reduces
to a half-sample replicate method.
}
\usage{
bootweights(strata, psu, replicates = 50, fpc = NULL,
         fpctype = c("population", "fraction", "correction"),
         compress = TRUE)
subbootweights(strata, psu, replicates = 50, compress = TRUE)
mrbweights(clusters, stratas, fpcs, replicates=50, 
         multicore=getOption("survey.multicore"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strata}{Identifier for sampling strata (top level only)}
  \item{stratas}{data frame of strata for all stages of sampling}
  \item{psu}{Identifier for primary sampling units}
  \item{clusters}{data frame of identifiers for sampling units at each stage}
  \item{replicates}{Number of bootstrap replicates}
  \item{fpc}{Finite population correction (top level only) }
  \item{fpctype}{Is \code{fpc} the population size, sampling fraction,
    or 1-sampling fraction?}
  \item{fpcs}{\code{survey_fpc} object with population and sample size at each stage}
  \item{compress}{Should the replicate weights be compressed?}
   \item{multicore}{Use the \code{multicore} package to generate the replicates in parallel}
}

\value{
 A set of replicate weights
}

\section{warning}{With \code{multicore=TRUE} the resampling procedure does not
use the current random seed, so the results cannot be exactly
reproduced even by using \code{set.seed()}}

\note{
These bootstraps are strictly appropriate only when the first stage of
sampling is a simple or stratified random sample of PSUs with or
without replacement, and not (eg) for PPS sampling.  The functions
will not enforce simple random sampling, so they can be used
(approximately) for data that have had non-response corrections and
other weight adjustments.  It is preferable to apply these adjustments
after creating the bootstrap replicate weights, but that may not be
possible with public-use data.

}
\references{Canty AJ, Davison AC. (1999) Resampling-based variance
  estimation for labour force surveys. The Statistician 48:379-391

  Judkins, D. (1990), "Fay's Method for Variance Estimation" Journal of Official Statistics, 6, 223-239.

  Preston J. (2009) Rescaled bootstrap for stratified multistage sampling. Survey Methodology 35(2) 227-234

  Rao JNK, Wu CFJ. Bootstrap inference for sample surveys. Proc Section
  on Survey Research Methodology. 1993 (866--871)
}

\seealso{\code{\link{as.svrepdesign}}}

\keyword{survey}% at least one, from doc/KEYWORDS

