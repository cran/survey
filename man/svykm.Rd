\name{svykm}
\alias{svykm}
\alias{plot.svykm}
\alias{plot.svykmlist}
\alias{lines.svykm}
\alias{quantile.svykm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate survival function. }
\description{
Estimates the survival function using a weighted Kaplan-Meier
estimator. Standard errors are not currently available.
}
\usage{
svykm(formula, design, ...)
\method{plot}{svykm}(x,xlab="time",ylab="Proportion surviving",
  ylim=c(0,1),...)
\method{lines}{svykm}(x,xlab="time",type="s",...)
\method{plot}{svykmlist}(x, pars=NULL, ...)
\method{quantile}{svykm}(x, probs=c(0.75,0.5,0.25),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Two-sided formula. The response variable should be a right-censored
    \code{Surv} object}
  \item{design}{survey design object}
  \item{\dots}{in \code{plot} and \code{lines} methods, graphical
parameters }
   \item{x}{a \code{svykm} or \code{svykmlist} object}
   \item{xlab,ylab,ylim,type}{as for \code{plot}}
   \item{pars}{A list of vectors of graphical parameters for the
separate curves in a \code{svykmlist} object}
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

biasmodel<-glm(I(trt>0)~age*edema,data=pbc)
pbc$randprob<-fitted(biasmodel)

dpbc<-svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,trt>0))

s1<-svykm(Surv(time,status)~1, design=dpbc)
s2<-svykm(Surv(time,status)~I(bili>6), design=dpbc)

plot(s1)
plot(s2)
plot(s2, lwd=2, pars=list(lty=c(1,2),col=c("purple","forestgreen")))

quantile(s1, probs=c(0.9,0.75,0.5,0.25,0.1))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{survival}% __ONLY ONE__ keyword per line
