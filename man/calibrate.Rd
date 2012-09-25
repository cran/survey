\name{calibrate}
\alias{calibrate.survey.design2}
\alias{calibrate.svyrep.design}
%\alias{is.calibrated}
\alias{calibrate}
%\alias{regcalibrate.survey.design2}
%\alias{regcalibrate.svyrep.design}
%\alias{regcalibrate}
\alias{calibrate.twophase}
\alias{grake}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calibration (GREG) estimators}
\description{
  Calibration, generalized raking, or GREG estimators generalise post-stratification and
  raking by calibrating a sample to the marginal totals of
  variables in a linear regression model.  This function reweights the
  survey design and adds additional information that is used by
  \code{svyrecvar} to reduce the estimated standard errors.
}
\usage{
calibrate(design,...)
\method{calibrate}{survey.design2}(design, formula, population,
       aggregate.stage=NULL, stage=0, variance=NULL,
       bounds=c(-Inf,Inf), calfun=c("linear","raking","logit"),
       maxit=50,epsilon=1e-7,verbose=FALSE,force=FALSE,trim=NULL,...)
\method{calibrate}{svyrep.design}(design, formula, population,compress=NA,
       aggregate.index=NULL, variance=NULL, bounds=c(-Inf,Inf),
       calfun=c("linear","raking","logit"),
       maxit=50, epsilon=1e-7, verbose=FALSE,force=FALSE,trim=NULL, ...)
\method{calibrate}{twophase}(design, phase=2,formula, population,
       calfun=c("linear","raking","logit","rrz"),...)
grake(mm,ww,calfun,eta=rep(0,NCOL(mm)),bounds,population,epsilon, verbose,maxit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{survey design object}
  \item{formula}{model formula for calibration model, or list of
    formulas for each margin}
  \item{population}{Vectors of population column totals for the model matrix in the
    calibration model, or list of such vectors for each
    cluster, or list of tables for each margin. Required except for two-phase designs}
  \item{compress}{compress the resulting replicate weights if
    \code{TRUE} or if \code{NA} and weights were previously compressed}
  \item{stage}{See Details below}
  \item{variance}{Coefficients for variance in calibration model (see
    Details below)}
  \item{aggregate.stage}{An integer. If not \code{NULL}, make calibration weights
    constant within sampling units at this stage.}
  \item{aggregate.index}{A vector or one-sided formula. If not \code{NULL}, make calibration weights
    constant within levels of this variable}
  \item{bounds}{Bounds for the calibration weights, optional
    except for \code{calfun="logit"}}
  \item{trim}{Weights outside this range will be trimmed to these bounds.}
  \item{...}{options for other methods}
  \item{calfun}{Calibration function: see below}
  \item{maxit}{Number of iterations}
  \item{epsilon}{tolerance in matching population total. Either a single
  number or a vector of the same length as \code{population}}
  \item{verbose}{print lots of uninteresting information}
  \item{force}{Return an answer even if the specified accuracy was not achieved}
  \item{phase}{Phase of a two-phase design to calibrate (only
    \code{phase=2} currently implemented.)}
  \item{mm}{model matrix}
  \item{ww}{vector of weights}
  \item{eta}{starting values for iteration}
}
\details{

  The \code{formula} argument specifies a model matrix, and the
  \code{population} argument is the population column sums of this
  matrix.

  For the important special case where the calibration totals are (possibly
  overlapping) marginal tables of factor variables, as in classical
  raking, the \code{formula} and \code{population} arguments may be
  lists in the same format as the input to \code{\link{rake}}.

  If the \code{population} argument has a names attribute it will be
  checked against the names produced by \code{model.matrix(formula)} and
  reordered if necessary.  This protects against situations where the
  (locale-dependent) ordering of factor levels is not what you expected.

  Numerical instabilities may result if the sampling weights in the
  \code{design} object are wrong by multiple orders of magnitude. The
  code now attempts to rescale the weights first, but it is better for
  the user to ensure that the scale is reasonable. 
  
  The \code{calibrate} function implements linear, bounded linear,
  raking, bounded raking, and logit calibration functions. All except
  unbounded linear calibration use the Newton-Raphson algorithm
  described by Deville et al (1993). This algorithm is exposed for other
  uses in the \code{grake} function.  Unbounded linear calibration uses
  an algorithm that is less sensitive to collinearity. The calibration
  function may be specified as a string naming one of the three built-in
  functions or as an object of class \code{calfun}, allowing
  user-defined functions. See \code{\link{make.calfun}} for details.

  Calibration with bounds, or on highly collinear data, may fail. If
  \code{force=TRUE} the approximately calibrated design object will
  still be returned (useful for examining why it failed). A failure in
  calibrating a set of replicate weights when the sampling weights were
  successfully calibrated will give only a warning, not an error.
  
  When calibration to the desired set of bounds is not possible, another option is
  to trim weights. To do this set \code{bounds} to a looser set of bounds
  for which calibration is achievable and set \code{trim} to the tighter
  bounds. Weights outside the bounds will be trimmed to the bounds, and
  the excess weight distributed over other observations in proportion to
  their sampling weight (and so this may put some other observations
  slightly over the trimming bounds). The projection matrix used in computing
  standard errors is based on the feasible bounds specified by the
  \code{bounds} argument.  See also \code{\link{trimWeights}},
  which trims the final weights in a design object rather than the
  calibration adjustments.

  
  For two-phase designs \code{calfun="rrz"} estimates the sampling
  probabilities using logistic regression as described by Robins et al
  (1994). \code{\link{estWeights}} will do the same thing.

  Calibration may result in observations within the last-stage sampling
  units having unequal weight even though they necessarily are sampled
  together.  Specifying \code{aggegrate.stage} ensures that the
  calibration weight adjustments are constant within sampling units at
  the specified stage; if the original sampling weights were equal the
  final weights will also be equal.  The algorithm is as described by
  Vanderhoeft (2001, section III.D). Specifying \code{aggregate.index}
  does the same thing for replicate weight designs; a warning will be
  given if the original weights are not constant within levels of
  \code{aggregate.index}.
  
  In a model with two-stage sampling, population totals may be available
  for the PSUs actually sampled, but not for the whole population.  In
  this situation, calibrating within each PSU reduces with second-stage
  contribution to variance. This generalizes to multistage sampling.
  The \code{stage} argument specifies which stage of sampling the totals
  refer to.  Stage 0 is full population totals, stage 1 is totals for
  PSUs, and so on.  The default, \code{stage=NULL} is interpreted as
  stage 0 when a single population vector is supplied and stage 1 when a
  list is supplied. Calibrating to PSU totals will fail (with a message
  about an exactly singular matrix) for PSUs that have fewer
  observations than the number of calibration variables.
  
  For unbounded linear calibration only, the variance in the calibration
  model may depend on covariates.  If \code{variance=NULL} the
  calibration model has constant variance.  If \code{variance} is not \code{NULL}
  it specifies a linear combination of the columns of the model matrix
  and the calibration variance is proportional to that linear
  combination.

  The design matrix specified by formula (after any aggregation) must be
  of full rank, with one exception. If the population total for a column
  is zero and all the observations are zero the column will be
  ignored. This allows the use of factors where the population happens
  to have no observations at some level.

  In a two-phase design, \code{population} may be omitted when
  \code{phase=2}, to specify calibration to the phase-one sample. If the
  two-phase design object was constructed using the more memory-efficient
  \code{method="approx"} argument to \code{\link{twophase}}, calibration of the first
  phase of sampling to the population is not supported.


}
\value{
  A survey design object.
}

\references{
Deville J-C, Sarndal C-E, Sautory O (1993) Generalized Raking
Procedures in Survey Sampling. JASA 88:1013-1020

Kalton G, Flores-Cervantes I (2003) "Weighting methods" J Official
 Stat 19(2) 81-97
  
Lumley T, Shaw PA, Dai JY (2011) "Connections between survey calibration estimators and semiparametric models for incomplete data" International Statistical Review. 79:200-220. (with discussion 79:221-232)

Sarndal C-E, Swensson B, Wretman J. "Model Assisted Survey
Sampling". Springer. 1991.

Rao JNK, Yung W, Hidiroglou MA (2002)   Estimating equations for the
analysis of survey data using poststratification information. Sankhya
64 Series A Part 2, 364-378.

Robins JM, Rotnitzky A, Zhao LP. (1994) Estimation of regression
coefficients when some regressors are not always observed. Journal of
the American Statistical Association, 89, 846-866.

Vanderhoeft C (2001) Generalized Calibration at Statistics
Belgium. Statistics Belgium Working Paper No 3.
\url{http://www.statbel.fgov.be/studies/paper03_en.asp}
}

\seealso{ \code{\link{postStratify}}, \code{\link{rake}} for other ways
  to use auxiliary information
  
 \code{\link{twophase}} and \code{vignette("epi")} for an example of calibration in two-phase designs

 \code{survey/tests/kalton.R} for examples replicating those in Kalton & Flores-Cervantes (2003)

 \code{\link{make.calfun}} for user-defined calibration distances.

 \code{\link{trimWeights}} to trim final weights rather than calibration adjustments.
}
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

pop.totals<-c(`(Intercept)`=6194, stypeH=755, stypeM=1018)

## For a single factor variable this is equivalent to
## postStratify

(dclus1g<-calibrate(dclus1, ~stype, pop.totals))

svymean(~api00, dclus1g)
svytotal(~enroll, dclus1g)
svytotal(~stype, dclus1g)

## Make weights constant within school district
(dclus1agg<-calibrate(dclus1, ~stype, pop.totals, aggregate=1))
svymean(~api00, dclus1agg)
svytotal(~enroll, dclus1agg)
svytotal(~stype, dclus1agg)


## Now add sch.wide
(dclus1g2 <- calibrate(dclus1, ~stype+sch.wide, c(pop.totals, sch.wideYes=5122)))

svymean(~api00, dclus1g2)
svytotal(~enroll, dclus1g2)
svytotal(~stype, dclus1g2)

## Finally, calibrate on 1999 API and school type

(dclus1g3 <- calibrate(dclus1, ~stype+api99, c(pop.totals, api99=3914069)))

svymean(~api00, dclus1g3)
svytotal(~enroll, dclus1g3)
svytotal(~stype, dclus1g3)


## Same syntax with replicate weights
rclus1<-as.svrepdesign(dclus1)

(rclus1g3 <- calibrate(rclus1, ~stype+api99, c(pop.totals, api99=3914069)))

svymean(~api00, rclus1g3)
svytotal(~enroll, rclus1g3)
svytotal(~stype, rclus1g3)

(rclus1agg3 <- calibrate(rclus1, ~stype+api99, c(pop.totals,api99=3914069), aggregate.index=~dnum))

svymean(~api00, rclus1agg3)
svytotal(~enroll, rclus1agg3)
svytotal(~stype, rclus1agg3)


###
## Bounded weights
range(weights(dclus1g3)/weights(dclus1))
dclus1g3b <- calibrate(dclus1, ~stype+api99, c(pop.totals, api99=3914069),bounds=c(0.6,1.6))
range(weights(dclus1g3b)/weights(dclus1))

svymean(~api00, dclus1g3b)
svytotal(~enroll, dclus1g3b)
svytotal(~stype, dclus1g3b)

## trimming
dclus1tr <- calibrate(dclus1, ~stype+api99, c(pop.totals, api99=3914069),bounds=c(0.5,2), trim=c(2/3,3/2))
svymean(~api00+api99+enroll, dclus1tr)
svytotal(~stype,dclus1tr)
range(weights(dclus1tr)/weights(dclus1))

rclus1tr <- calibrate(rclus1, ~stype+api99, c(pop.totals, api99=3914069),bounds=c(0.5,2), trim=c(2/3,3/2))
svymean(~api00+api99+enroll, rclus1tr)
svytotal(~stype,rclus1tr)

## Input in the same format as rake() for classical raking
pop.table <- xtabs(~stype+sch.wide,apipop)
pop.table2 <- xtabs(~stype+comp.imp,apipop)
dclus1r<-rake(dclus1, list(~stype+sch.wide, ~stype+comp.imp),
               list(pop.table, pop.table2))
gclus1r<-calibrate(dclus1, formula=list(~stype+sch.wide, ~stype+comp.imp), population=list(pop.table, pop.table2),calfun="raking")
svymean(~api00+stype, dclus1r)
svymean(~api00+stype, gclus1r)


## generalised raking
dclus1g3c <- calibrate(dclus1, ~stype+api99, c(pop.totals,
    api99=3914069), calfun="raking")
range(weights(dclus1g3c)/weights(dclus1))

(dclus1g3d <- calibrate(dclus1, ~stype+api99, c(pop.totals,
    api99=3914069), calfun=cal.logit, bounds=c(0.5,2.5)))
range(weights(dclus1g3d)/weights(dclus1))



## Ratio estimators are calibration estimators
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
svytotal(~api.stu,dstrat)

common<-svyratio(~api.stu, ~enroll, dstrat, separate=FALSE)
predict(common, total=3811472)

pop<-3811472
## equivalent to (common) ratio estimator
dstratg1<-calibrate(dstrat,~enroll-1, pop, variance=1)
svytotal(~api.stu, dstratg1)


}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
