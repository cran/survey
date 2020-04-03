\name{svymle}
\alias{svymle}
\alias{print.svymle}
\alias{coef.svymle}
\alias{summary.svymle}
\alias{vcov.svymle}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Maximum pseudolikelihood estimation in complex surveys}
\description{
Maximises a user-specified likelihood parametrised by multiple linear
predictors to data from a complex sample survey and computes the
sandwich variance estimator of the coefficients. Note that this function
maximises an estimated population likelihood, it is not the sample MLE.
}
\usage{
svymle(loglike, gradient = NULL, design, formulas, start = NULL, control
= list(), na.action="na.fail", method=NULL, lower=NULL,upper=NULL,influence=FALSE,...)
\method{summary}{svymle}(object, stderr=c("robust", "model"),...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{loglike}{vectorised loglikelihood function}
  \item{gradient}{Derivative of \code{loglike}. Required for variance computation and helpful for fitting}
  \item{design}{ a \code{survey.design} object }
  \item{formulas}{A list of formulas specifying the variable and linear predictors: see Details below}
  \item{start}{Starting values for parameters}
  \item{control}{control options for the optimiser: see the help page
    for the optimiser you are using.}
  \item{lower,upper}{Parameter bounds for \code{bobyqa}}
  \item{influence}{Return the influence functions (primarily for svyby)}
  \item{na.action}{Handling of \code{NA}s}
  \item{method}{\code{"nlm"} to use \code{nlm},\code{"uobyqa"} or
    \code{"bobyqa"} to use those optimisers from the \code{minqa}
    package;  otherwise passed to \code{\link{optim}}}
  \item{\dots}{Arguments to \code{loglike} and \code{gradient} that are
    not to be optimised over.}
  \item{object}{\code{svymle} object}
  \item{stderr}{Choice of standard error estimator. The default is a
    standard sandwich estimator. See Details below.}
}
\details{
Optimization is done by \code{\link{nlm}} by default or if
\code{method=="nlm"}. Otherwise \code{\link{optim}} is used and \code{method}
specifies the method and \code{control} specifies control parameters.

  The \code{design} object contains all the data and design information
from the survey, so all the formulas refer to variables in this object.
The \code{formulas} argument needs to specify the response variable and
a linear predictor for each freely varying argument of \code{loglike}.

Consider for example the \code{\link{dnorm}} function, with arguments
\code{x}, \code{mean}, \code{sd} and \code{log}, and suppose we want to
estimate the mean of \code{y} as a linear function of a variable
\code{z}, and to estimate a constant standard deviation.  The \code{log}
argument must be fixed at \code{FALSE} to get the loglikelihood.  A
\code{formulas} argument would be \code{list(~y, mean=~z, sd=~1)}. Note
that the data variable \code{y} must be the first argument to
\code{dnorm} and the first formula and that all the other formulas are
labelled.  It is also permitted to have the data variable as the
left-hand side of one of the formulas: eg \code{list( mean=y~z, sd=~1)}.

The two optimisers from the \code{minqa} package do not use any
derivatives to be specified for optimisation, but they do assume
that the function is smooth enough for a quadratic approximation, ie,
that two derivatives exist.


The usual variance estimator for MLEs in a survey sample is a `sandwich'
variance that requires the score vector and the information matrix. It
requires only sampling assumptions to be valid (though some model
assumptions are required for it to be useful). This is the
\code{stderr="robust"} option, which is available only when the \code{gradient}
argument was specified.

If the model is correctly specified and the sampling is at random
conditional on variables in the model then standard errors based on just
the information matrix will be approximately valid.  In particular, for
independent sampling where weights and strata depend on variables in the
model the \code{stderr="model"} should work fairly well.
}
\value{
  An object of class \code{svymle}
}
\author{Thomas Lumley}

\seealso{\code{\link{svydesign}}, \code{\link{svyglm}}}

\examples{

 data(api)

 dstrat<-svydesign(id=~1, strata=~stype, weight=~pw, fpc=~fpc, data=apistrat)

 ## fit with glm
 m0 <- svyglm(api00~api99+ell,family="gaussian",design=dstrat)
 ## fit as mle (without gradient)
 m1 <- svymle(loglike=dnorm,gradient=NULL, design=dstrat, 
    formulas=list(mean=api00~api99+ell, sd=~1),
    start=list(c(80,1,0),c(20)), log=TRUE)
 ## with gradient
 gr<- function(x,mean,sd,log){
	 dm<-2*(x - mean)/(2*sd^2)
	 ds<-(x-mean)^2*(2*(2 * sd))/(2*sd^2)^2 - sqrt(2*pi)/(sd*sqrt(2*pi))
         cbind(dm,ds)
      }
 m2 <- svymle(loglike=dnorm,gradient=gr, design=dstrat, 
    formulas=list(mean=api00~api99+ell, sd=~1),
    start=list(c(80,1,0),c(20)), log=TRUE, method="BFGS")

 summary(m0)
 summary(m1,stderr="model")
 summary(m2)

 ## Using offsets
 m3 <- svymle(loglike=dnorm,gradient=gr, design=dstrat, 
    formulas=list(mean=api00~api99+offset(ell)+ell, sd=~1),
    start=list(c(80,1,0),c(20)), log=TRUE, method="BFGS")


## demonstrating multiple linear predictors

 m3 <- svymle(loglike=dnorm,gradient=gr, design=dstrat, 
    formulas=list(mean=api00~api99+offset(ell)+ell, sd=~stype),
    start=list(c(80,1,0),c(20,0,0)), log=TRUE, method="BFGS")



 ## More complicated censored lognormal data example
 ## showing that the response variable can be multivariate

 data(pbc, package="survival")
 pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
 biasmodel<-glm(randomized~age*edema,data=pbc)
 pbc$randprob<-fitted(biasmodel)
 dpbc<-svydesign(id=~1, prob=~randprob, strata=~edema,
    data=subset(pbc,randomized))


## censored logNormal likelihood
 lcens<-function(x,mean,sd){
    ifelse(x[,2]==1,
           dnorm(log(x[,1]),mean,sd,log=TRUE),
           pnorm(log(x[,1]),mean,sd,log=TRUE,lower.tail=FALSE)
           )
 }

 gcens<- function(x,mean,sd){

        dz<- -dnorm(log(x[,1]),mean,sd)/pnorm(log(x[,1]),mean,sd,lower.tail=FALSE)

        dm<-ifelse(x[,2]==1,
                   2*(log(x[,1]) - mean)/(2*sd^2),
                   dz*-1/sd)
        ds<-ifelse(x[,2]==1,
                   (log(x[,1])-mean)^2*(2*(2 * sd))/(2*sd^2)^2 - sqrt(2*pi)/(sd*sqrt(2*pi)),
                   ds<- dz*-(log(x[,1])-mean)/(sd*sd))
        cbind(dm,ds)      
 }

m<-svymle(loglike=lcens, gradient=gcens, design=dpbc, method="newuoa",
      formulas=list(mean=I(cbind(time,status>0))~bili+protime+albumin,
                    sd=~1),
         start=list(c(10,0,0,0),c(1)))

summary(m)

## the same model, but now specifying the lower bound of zero on the
## log standard deviation

mbox<-svymle(loglike=lcens, gradient=gcens, design=dpbc, method="bobyqa",
      formulas=list(mean=I(cbind(time,status>0))~bili+protime+albumin, sd=~1),
      lower=list(c(-Inf,-Inf,-Inf,-Inf),0), upper=Inf,
      start=list(c(10,0,0,0),c(1)))


## The censored lognormal model is now available in svysurvreg()

summary(svysurvreg(Surv(time,status>0)~bili+protime+albumin,
        design=dpbc,dist="lognormal"))

## compare svymle scale value after log transformation
svycontrast(m, quote(log(`sd.(Intercept)`)))



}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{optimize}% __ONLY ONE__ keyword per line
