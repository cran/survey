\name{estweights}
\alias{estWeights}
\alias{estWeights.twophase}
\alias{estWeights.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimated weights for missing data}
\description{
Creates or adjusts a two-phase survey design object using a logistic
regression model for second-phase sampling probability.  This function
should be particularly useful in reweighting to account for missing data.  }
\usage{
estWeights(data,formula,...)
\method{estWeights}{twophase}(data,formula,...)
\method{estWeights}{data.frame}(data,formula, subset=NULL, strata=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{twophase design object or data frame}
  \item{formula}{Predictors for estimating weights}
  \item{subset}{Subset of data frame with complete data (ie phase 1).
    If \code{NULL} use all complete cases}
  \item{strata}{Stratification (if any) of phase 2 sampling}
  \item{\dots}{for future expansion}
}
\details{
  If \code{data} is a data frame, \code{estWeights} first creates a
  two-phase design object. The \code{strata} argument is used only to
  compute finite population corrections, the same variables must be
  included in \code{formula} to compute stratified sampling probabilities.
 
  With a two-phase design object, \code{estWeights} estimates the sampling
  probabilities using logistic regression as described by Robins et al
  (1994) and adds information to the object to enable correct sandwich
  standard errors to be computed.
  
  The effect on a two-phase design object is very similar to
  \code{\link{calibrate}}, and is identical when \code{formula}
  specifies a saturated model.
}
\value{
  A two-phase survey design object.
}

\references{

Robins JM, Rotnitzky A, Zhao LP. (1994) Estimation of regression
coefficients when some regressors are not always observed. Journal of
the American Statistical Association, 89, 846-866.

}

\seealso{ \code{\link{postStratify}},
  \code{\link{calibrate}}, \code{\link{twophase}}}
\examples{
data(airquality)

## ignoring missingness, using model-based standard error
summary(lm(log(Ozone)~Temp+Wind, data=airquality))

## Without covariates to predict missingness we get
## same point estimates, but different (sandwich) standard errors
daq<-estWeights(airquality, formula=~1,subset=~I(!is.na(Ozone)))
summary(svyglm(log(Ozone)~Temp+Wind,design=daq))

## Reweighting based on weather, month
d2aq<-estWeights(airquality, formula=~Temp+Wind+Month,
                 subset=~I(!is.na(Ozone)))
summary(svyglm(log(Ozone)~Temp+Wind,design=d2aq))

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
