\name{svykm}
\alias{svykm}
\alias{plot.svykm}
\alias{plot.svykmlist}
\alias{lines.svykm}
\alias{quantile.svykm}
\alias{confint.svykm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate survival function. }
\description{
Estimates the survival function using a weighted Kaplan-Meier
estimator. 
}
\usage{
svykm(formula, design,se=FALSE, ...)
\method{plot}{svykm}(x,xlab="time",ylab="Proportion surviving",
  ylim=c(0,1),ci=NULL,lty=1,...)
\method{lines}{svykm}(x,xlab="time",type="s",ci=FALSE,lty=1,...)
\method{plot}{svykmlist}(x, pars=NULL, ci=FALSE,...)
\method{quantile}{svykm}(x, probs=c(0.75,0.5,0.25),ci=FALSE,level=0.95,...)
\method{confint}{svykm}(object,parm,level=0.95,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Two-sided formula. The response variable should be a right-censored
    \code{Surv} object}
  \item{design}{survey design object}
  \item{se}{Compute standard errors? This is slow for moderate to large
    data sets}
  \item{\dots}{in \code{plot} and \code{lines} methods, graphical
parameters }
   \item{x}{a \code{svykm} or \code{svykmlist} object}
   \item{xlab,ylab,ylim,type}{as for \code{plot}}
   \item{lty}{Line type, see \code{\link{par}}}
   \item{ci}{Plot (or return, for\code{quantile}) the confidence interval}
   \item{pars}{A list of vectors of graphical parameters for the
     separate curves in a \code{svykmlist} object}
   \item{object}{A \code{svykm} object}
   \item{parm}{vector of times to report confidence intervals}
   \item{level}{confidence level}
   \item{probs}{survival probabilities for computing survival quantiles
     (note that these are the complement of the usual
     \code{\link{quantile}} input, so 0.9 means 90\% surviving, not 90\% dead)}
}
\value{
  For \code{svykm}, an object of class \code{svykm} for a single curve or \code{svykmlist}
  for multiple curves.
}
\details{
  When standard errors are computed, the survival curve is
  actually the Aalen (hazard-based) estimator rather than the
  Kaplan-Meier estimator. 
  
  The standard error computations use memory proportional to the sample
  size times the square of the number of events. This can be a lot.
  
  In the case of equal-probability cluster sampling without replacement
  the computations are essentially the same as those of Williams (1995),
  and the same linearization strategy is used for other designs.

  Confidence intervals are computed on the log(survival) scale,
  following the default in \code{survival} package, which was based on
  simulations by Link(1984).

  Confidence intervals for quantiles use Woodruff's method: the interval
  is the intersection of the horizontal line at the specified quantile
  with the pointwise confidence band around the survival curve.
}
\references{
Link, C. L. (1984). Confidence intervals for the survival function using
Cox's proportional hazards model with covariates. Biometrics 40,
601-610.

Williams RL (1995) "Product-Limit Survival Functions with Correlated 
Survival Times" Lifetime Data Analysis 1: 171--186

Woodruff RS (1952) Confidence intervals for medians and other
position measures. JASA 57, 622-627.  
 }
 \seealso{\code{\link{predict.svycoxph}} for survival curves from a Cox model
 }
\examples{
data(pbc, package="survival")
pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
biasmodel<-glm(randomized~age*edema,data=pbc)
pbc$randprob<-fitted(biasmodel)

dpbc<-svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,randomized))

s1<-svykm(Surv(time,status>0)~1, design=dpbc)
s2<-svykm(Surv(time,status>0)~I(bili>6), design=dpbc)

plot(s1)
plot(s2)
plot(s2, lwd=2, pars=list(lty=c(1,2),col=c("purple","forestgreen")))

quantile(s1, probs=c(0.9,0.75,0.5,0.25,0.1))

s3<-svykm(Surv(time,status>0)~I(bili>6), design=dpbc,se=TRUE)
plot(s3[[2]],col="purple")

confint(s3[[2]], parm=365*(1:5))
quantile(s3[[1]], ci=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{survival}% __ONLY ONE__ keyword per line
