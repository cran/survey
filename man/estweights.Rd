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
\method{estWeights}{twophase}(data,formula=NULL, working.model=NULL,...)
\method{estWeights}{data.frame}(data,formula=NULL, working.model=NULL,
      subset=NULL, strata=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{twophase design object or data frame}
  \item{formula}{Predictors for estimating weights}
  \item{working.model}{Model fitted to complete (ie phase 1) data}
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

  An alternative to specifying \code{formula} is to specify
  \code{working.model}. The estimating functions from this model will be
  used as predictors of the sampling probabilities, which will increase
  efficiency to the extent that the working model and the model of
  interest estimate the same parameters (Kulich \& Lin 2004).
  
  The effect on a two-phase design object is very similar to
  \code{\link{calibrate}}, and is identical when \code{formula}
  specifies a saturated model.
}
\value{
  A two-phase survey design object.
}

\references{
Breslow NE, Lumley T, Ballantyne CM, Chambless LE, Kulich M. (2009) Using the Whole Cohort in the Analysis of Case-Cohort Data.  Am J Epidemiol. 2009 Jun 1;169(11):1398-405.

Robins JM, Rotnitzky A, Zhao LP. (1994) Estimation of regression
coefficients when some regressors are not always observed. Journal of
the American Statistical Association, 89, 846-866.

Kulich M, Lin DY (2004). Improving the Efficiency of Relative-Risk
Estimation in Case-Cohort Studies. Journal of the American Statistical Association, Vol. 99,  pp.832-844 

Lumley T, Shaw PA, Dai JY (2011) "Connections between survey calibration estimators and semiparametric models for incomplete data" International Statistical Review. 79:200-220. (with discussion 79:221-232)

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
