\name{calibrate}
\alias{calibrate.survey.design2}
\alias{calibrate.svyrep.design}
\alias{is.calibrated}
\alias{calibrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{G-calibration (GREG) estimators}
\description{
  G-calibration (GREG) estimators generalise post-stratification and
  raking by calibrating a sample to the marginal totals of
  variables in a linear regression model.  This function reweights the
  survey design and adds additional information that is used by
  \code{svyrecvar} to reduce the estimated standard errors.
}
\usage{
calibrate(design,...)
\method{calibrate}{survey.design2}(design, formula, population,stage=NULL,lambda=NULL,...)
\method{calibrate}{svyrep.design}(design, formula, population, compress=NA,lambda=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{survey design object}
  \item{formula}{model formula for calibration model}
  \item{population}{Vectors of population column totals for the model matrix in the
    calibration model, or list of such vectors for each cluster.}
  \item{compress}{compress the resulting replicate weights if
    \code{TRUE} or if \code{NA} and weights were previously compressed}
  \item{stage}{See Details below}
  \item{lambda}{Coefficients for variance in calibration model (see
    Details below)}
  \item{...}{options for other methods}
}
\details{
  In a model with two-stage sampling population totals may be available
  for the PSUs actually sampled, but not for the whole population.  In
  this situation, calibrating within each PSU reduces with second-stage
  contribution to variance. This generalizes to multistage sampling.
  
  The \code{stage} argument specifies which stage of sampling the
  totals refer to.  Stage 0 is full population totals, stage 1 is
  totals for PSUs, and so on.  The default, \code{stage=NULL} is
  interpreted as stage 0 when a single population vector is supplied
  and stage 1 when a list is supplied.
  
  If \code{lambda=NULL} the calibration model has constant
  variance. The model must explicitly or implicitly contain an
  intercept. If \code{lambda} is not \code{NULL} it specifies a linear
  combination of the columns of the model matrix and the calibration
  variance is proportional to that linear combination.
}
\value{
  A survey design object.
}

\references{
Sarndal CA, Swensson B, Wretman J. "Model Assisted Survey
Sampling". Springer. 1991.


Rao JNK, Yung W, Hidiroglou MA (2002)   Estimating equations for the
analysis of survey data using poststratification information. Sankhya
64 Series A Part 2, 364-378.
}

\seealso{ \code{\link{postStratify}}, \code{\link{rake}}}
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

## Ratio estimators
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
rstrat<-as.svrepdesign(dstrat)

svytotal(~api.stu,dstrat)

common<-svyratio(~api.stu, ~enroll, dstrat, separate=FALSE)
sep<-svyratio(~api.stu,~enroll, dstrat,separate=TRUE)
stratum.totals<-list(E=1877350, H=1013824, M=920298)
predict(sep, total=stratum.totals)
predict(common, total=do.call("sum",stratum.totals))

pop<-colSums(model.matrix(~stype*enroll-1,model.frame(~stype*enroll,apipop)))
pop
## common ratio
dstratg1<-calibrate(dstrat,~enroll-1, pop[4], lambda=1)
svytotal(~api.stu, dstratg1)
rstratg1<-calibrate(rstrat,~enroll-1, pop[4], lambda=1)
svytotal(~api.stu, rstratg1)

## similar (but not identical) to separate ratio.
dstratg2<-calibrate(dstrat,~stype*enroll-1, pop,lambda=c(0,0,0,1,0,0))
svytotal(~api.stu,dstratg2)

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
