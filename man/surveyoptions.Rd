\name{surveyoptions}
\alias{survey.lonely.psu}
\alias{survey.ultimate.cluster}
\alias{survey.adjust.domain.lonely}
\alias{survey.want.obsolete}
\title{Options for the survey package}
\description{
  This help page documents the options that control the behaviour of
  the survey package.
  }
  
  \details{
    All the options for the survey package have names beginning with
    "survey".  Three of them control standard error estimation.
    
    When \code{options("survey.ultimate.cluster")} is \code{TRUE},
    standard error estimation is based on independence of PSUs at the
    first stage of sampling, without using any information about
    subsequent stages (this is what Stata and many other packages
    do). When \code{FALSE}, finite population corrections and variances
    are estimated recursively. See \code{\link{svyrecvar}} for more
    information.
    
    Handling of strata with a single PSU that are not certainty PSUs is
    controlled by \code{options("survey.lonely.psu")}. The default
    setting is \code{"fail"}, which gives an error. Use \code{"remove"}
    to ignore that PSU for variance computation, \code{"adjust"} to
    center the stratum at the population mean rather than the stratum
    mean, and \code{"average"} to replace the variance contribution of
    the stratum by the average variance contribution across strata.
    
    The variance formulas for domain estimation give well-defined,
    positive results when a stratum contains only one PSU with
    observations in the domain, but are not unbiased.  If
    \code{options("survey.adjust.domain.lonely")} is \code{TRUE} and
    \code{options("survey.lonely.psu")} is \code{"average"} or
    \code{"adjust"} the same adjustment for lonely PSUs will be used
    within a domain.
    
    The fourth option is \code{options("survey.want.obsolete")}. This
    controls the warnings about using the deprecated pre-2.9.0 survey
    design objects.
    
  }
  
  \keyword{survey}