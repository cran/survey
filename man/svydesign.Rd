\name{svydesign}
\alias{svydesign}
\alias{oldsvydesign}
\alias{summary.survey.design}
\alias{summary.survey.design2}
\alias{print.summary.survey.design}
\alias{print.summary.survey.design2}
\alias{print.survey.design}
\alias{print.survey.design2}
\alias{[<-.survey.design}
\alias{[.survey.design2}
\alias{model.frame.survey.design}
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
  \item{nest}{If \code{TRUE}, relabel cluster ids to enforce nesting
    within strata}
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
   population size is the sampling units.  That is, sampling 100 units
   from a population stratum of size 500 can be specified as 100 or as
   100/500=0.2.
   
   If population sizes are specified but not sampling probabilities or
   weights, the sampling probabilities will be computed from the
   population sizes assuming simple random sampling within strata.
   
   For multistage sampling the \code{id} argument should specify a
   formula with the cluster identifiers at each stage.  If subsequent
   stages are stratified \code{strata} should also be specified as a
   formula with stratum identifiers at each stage.  The population size
   for each level of sampling should also be specified in \code{fpc}.
   If \code{fpc} is not specified then sampling is assumed to be with
   replacement at the top level and only the first stage of cluster is
   used in computing variances. If \code{fpc} is specified but for fewer
   stages than \code{id}, sampling is assumed to be complete for
   subsequent stages.   The variance calculations for
   multistage sampling assume simple or stratified random sampling
   within clusters at each stage except possibly the last.
   
   
The \code{dim}, \code{"["}, \code{"[<-"} and na.action methods for
\code{survey.design} objects operate on the dataframe specified by
\code{variables} and ensure that the design information is properly
updated to correspond to the new data frame.  With the \code{"[<-"}
method the new value can be a \code{survey.design} object instead of a
data frame, but only the data frame is used. See also
\code{\link{subset.survey.design}} for a simple way to select
subpopulations.

The \code{model.frame} method extracts the observed data.


 If the strata with one only PSU are not self-representing (or they are,
but \code{svydesign} cannot tell based on \code{fpc}) then the handling
of these strata for variance computation is determined by
\code{options("survey.lonely.psu")}.  See \code{\link{svyCprod}} for
details.


}

\note{
  Use \code{oldsvydesign}, which has the same arguments, to create
  objects with the structure used before version 2.9.
  }

\value{
An object of class \code{survey.design}.
}
\author{Thomas Lumley}


\seealso{
  \code{\link{postStratify}} for post-stratification,
  \code{\link{as.svrepdesign}} for converting to replicate weight designs,
  \code{\link{subset.survey.design}} for domain estimates,
  \code{\link{update.survey.design}} to add variables.
}



\examples{
  data(api)
# stratified sample
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
# one-stage cluster sample
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
# two-stage cluster sample: weights computed from population sizes.
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)

## multistage sampling has no effect when fpc is not given, so
## these are equivalent.
dclus2wr<-svydesign(id=~dnum+snum, weights=weights(dclus2), data=apiclus2)
dclus2wr2<-svydesign(id=~dnum, weights=weights(dclus2), data=apiclus2)


## syntax for stratified cluster sample
##(though the data weren't really sampled this way)
svydesign(id=~dnum, strata=~stype, weights=~pw, data=apistrat, nest=TRUE)

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
\keyword{manip}% __ONLY ONE__ keyword per line
