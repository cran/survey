\name{svydesign}
\alias{svydesign}
\alias{summary.survey.design}
\alias{print.summary.survey.design}
\alias{print.survey.design}
\alias{[.survey.design}
\alias{[<-.survey.design}
\alias{na.omit.survey.design}
\alias{na.exclude.survey.design}
\alias{na.fail.survey.design}
\alias{dim.survey.design}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Survey sample analysis.}
\description{
  Specify a complex survey design.
}
\usage{
svydesign(ids, probs=NULL, strata = NULL, variables = NULL, fpc=NULL,
data = NULL, nest = FALSE, check.strata = !nest, weights=NULL) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{ids}{Formula or data frame specifying cluster ids from largest
    level to smallest level, \code{~0} or \code{~1} is a formula for no clusters.}
  \item{probs}{Formula or data frame specifying cluster sampling probabilities}
  \item{strata}{Formula or vector specifying strata, use \code{NULL} for no strata}
  \item{variables}{Formula or data frame specifying the variables
    measured in the survey. If \code{NULL}, the \code{data} argument is
    used.}
  \item{fpc}{Finite population correction: see Details below}
  \item{weights}{Formula or vector specifying sampling weights as an
    alternative to \code{prob}}
  \item{data}{Data frame to look up variables in the formula arguments}
  \item{nest}{If \code{TRUE}, relabel cluster ids to enforce nesting, eg
  if ids at second level of sampling are reused within first-level
  units}
  \item{check.strata}{If \code{TRUE}, check that clusters are nested in strata}.
}
\details{
   When analysing data from a complex survey, observations must be
   weighted inversely to their sampling probabilities, and the effects
   of stratification and of correlation induced by cluster sampling must
   be incorporated in standard errors.

   The \code{svydesign} object combines a data frame and all the survey
   design information needed to analyse it.  These objects are used by
   the survey modelling and summary functions.

   The finite population correction is used to reduce the variance when
   a substantial fraction of the total population of interest has been
   sampled. It may not be appropriate if the target of inference is the
   process generating the data rather than the statistics of a
   particular finite population.

   The finite population correction can be specified either as the total
   population size in each stratum or as the fraction of the total
   population that has been sampled. In either case the relevant
   population size is `primary sampling units', the largest clusters.
   That is, sampling 100 units from a population stratum of size 500 can
   be specified as 100 or as 100/500=0.2.  The finite population
   correction can be specified by a vector with one element for each
   individual (in which case it is an error for it to vary within a
   stratum) or as a data frame with one row per stratum.  The first
   column of the data frame should be a factor with the same levels as
   \code{strata} and the second column the finite population correction.
   
   If population sizes are specified but not sampling probabilities or
   weights, the sampling probabilities will be computed from the
   population sizes assuming simple random sampling within strata.

The \code{dim}, \code{"["}, \code{"[<-"} and na.action methods for
\code{survey.design} objects operate on the dataframe specified by
\code{variables} and ensure that the design information is properly
updated to correspond to the new data frame.  With the \code{"[<-"}
method the new value can be a \code{survey.design} object instead of a
data frame, but only the data frame is used. See also
\code{\link{subset.survey.design}} for a simple way to select
subpopulations.

The value of \code{options("survey.lonely.psu")} controls what happens
to strata containing only one cluster (PSU).See \code{\link{svyCprod}}
for details, especially if you have self-representing ("certainty") PSUs.

}
\value{
An object of class \code{survey.design}.
}
\author{Thomas Lumley}


\seealso{ \code{\link{svyglm}}, \code{\link{svymean}}, \code{\link{svyvar}}, \code{\link{svytable}}, \code{\link{svyquantile}},
\code{\link{subset.survey.design},
\code{\link{update.survey.design}}}
}

\examples{
  data(api)
# stratified sample
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
# one-stage cluster sample
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
# two-stage cluster sample
dclus2<-svydesign(id=~dnum+snum, weights=~pw, data=apiclus2)

## syntax for stratified cluster sample
##(though the data weren't really sampled this way)
svydesign(id=~dnum, strata=~stype, weights=~pw, data=apistrat, nest=TRUE)

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
\keyword{manip}% __ONLY ONE__ keyword per line
