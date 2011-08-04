\name{surveyoptions}
\alias{surveyoptions}
\alias{survey.lonely.psu}
\alias{survey.ultimate.cluster}
\alias{survey.adjust.domain.lonely}
\alias{survey.want.obsolete}
\alias{survey.drop.replicates}
\alias{survey.multicore}
\alias{survey.replicates.mse}
\title{Options for the survey package}
\description{
  This help page documents the options that control the behaviour of
  the survey package.
  }
  
  \details{
    All the options for the survey package have names beginning with
    "survey".  Four of them control standard error estimation.

    \code{options("survey.replicates.mse")} controls the default in
    \code{svrepdesign} and \code{as.svrepdesign} for computing
    variances. When \code{options("survey.replicates.mse")} is
    \code{TRUE}, the default is to create replicate weight designs that
    compute variances centered at the point estimate, rather than at the
    mean of the replicates.  The option can be overridden by specifying
    the \code{mse} argument explicitly in \code{svrepdesign} and
    \code{as.svrepdesign}. The default is \code{FALSE}.
  
    When \code{options("survey.ultimate.cluster")} is \code{TRUE},
    standard error estimation is based on independence of PSUs at the
    first stage of sampling, without using any information about
    subsequent stages. When \code{FALSE}, finite population corrections
    and variances are estimated recursively. See \code{\link{svyrecvar}}
    for more information. This option makes no difference unless
    first-stage finite population corrections are specified, in which
    case setting the option to \code{TRUE} gives the wrong answer for a
    multistage study. The only reason to use \code{TRUE} is for
    compatibility with other software that gives the wrong answer.
    
    Handling of strata with a single PSU that are not certainty PSUs is
    controlled by \code{options("survey.lonely.psu")}. The default
    setting is \code{"fail"}, which gives an error. Use \code{"remove"}
    to ignore that PSU for variance computation, \code{"adjust"} to
    center the stratum at the population mean rather than the stratum
    mean, and \code{"average"} to replace the variance contribution of
    the stratum by the average variance contribution across strata. As
    of version 3.4-2 \code{as.svrepdesign} also uses this option.
    
    The variance formulas for domain estimation give well-defined,
    positive results when a stratum contains only one PSU with
    observations in the domain, but are not unbiased.  If
    \code{options("survey.adjust.domain.lonely")} is \code{TRUE} and
    \code{options("survey.lonely.psu")} is \code{"average"} or
    \code{"adjust"} the same adjustment for lonely PSUs will be used
    within a domain. Note that this adjustment is not available for
    replicate-weight designs, nor (currently) for raked,
    post-stratified, or calibrated designs.
    
    The fourth option is \code{options("survey.want.obsolete")}. This
    controls the warnings about using the deprecated pre-2.9.0 survey
    design objects.

    The behaviour of replicate-weight designs for self-representing
    strata is controlled by \code{options("survey.drop.replicates")}.
    When \code{TRUE}, various optimizations are used that take advantage
    of the fact that these strata do not contribute to the variance.
    The only reason ever to use \code{FALSE} is if there is a bug in
    the code for these optimizations.

    The fifth option controls the use of multiple processors with the
    \code{multicore} package. This option should not affect the values
    computed by any of the survey functions.  If \code{TRUE}, all
    functions that are able to use multiple processors will do so by
    default. Using multiple processors may speed up calculations, but
    need not, especially if the computer is short on memory. The best
    strategy is probably to experiment with explicitly requesting
    \code{multicore=TRUE} in functions that support it, to see if there
    is an increase in speed before setting the global option. 
    
  }
  
  \keyword{survey}
