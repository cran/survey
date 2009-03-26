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
estimator. Standard errors are not currently available.
}
\usage{
svykm(formula, design,se=FALSE, ...)
\method{plot}{svykm}(x,xlab="time",ylab="Proportion surviving",
  ylim=c(0,1),ci=NULL,lty=1,...)
\method{lines}{svykm}(x,xlab="time",type="s",ci=FALSE,lty=1,...)
\method{plot}{svykmlist}(x, pars=NULL, ci=FALSE,...)
\method{quantile}{svykm}(x, probs=c(0.75,0.5,0.25),...)
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
   \item{ci}{Plot the confidence interval}
   \item{pars}{A list of vectors of graphical parameters for the
     separate curves in a \code{svykmlist} object}
   \item{object}{A \code{svykm} object}
   \item{parm}{vector of times to report confidence intervals}
   \item{level}{confidence level}
   \item{probs}{survival probabilities for computing survival quantiles
     (note that these are the complement of the usual \code{\link{quantile}} input}
}
\value{
  For \code{svykm}, an object of class \code{svykm} for a single curve or \code{svykmlist}
  for multiple curves.
}

\seealso{\code{\link{svycoxph}} }
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

confint(s3[[2]], 365*(1:5))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{survival}% __ONLY ONE__ keyword per line
