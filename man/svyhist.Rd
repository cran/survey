\name{svyhist}
\alias{svyhist}
\alias{svyboxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Histograms and boxplots}
\description{
Histograms and boxplots weighted by the sampling weights. 
}
\usage{
svyhist(formula, design, breaks = "Sturges",
      include.lowest = TRUE, right = TRUE, xlab = NULL,
       main = NULL, probability = TRUE, freq = !probability, ...)
svyboxplot(formula, design, all.outliers=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{One-sided formula for \code{svyhist}, two-sided for \code{svyboxplot}}
  \item{design}{A survey design object}
  \item{xlab}{x-axis label}
  \item{main}{Main title}
  \item{probability,freq}{Y-axis is probability density or frequency}
  \item{all.outliers}{Show all outliers in the boxplot, not just extremes}
  \item{breaks, include.lowest, right}{As for \code{\link{hist}}}
  \item{\dots}{Other arguments to \code{\link{hist}} or \code{\link{bxp}}}
}

\details{
The histogram breakpoints are computed as if the sample were a
simple random sample of the same size.

The grouping variable in \code{svyboxplot}, if present, must be a factor.

The boxplot whiskers go to the maximum and minimum observations or to
1.5 interquartile ranges beyond the end of the box, whichever is
closer. The maximum and minimum are plotted as outliers if they are
beyond the ends of the whiskers, but other outlying points are not
plotted unless \code{all.outliers=TRUE}.  \code{svyboxplot}
requires a two-sided formula; use \code{variable~1} for a single boxplot.
 }

\seealso{ \code{\link{svyplot}}}
\examples{
data(api)
dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat, 
    fpc = ~fpc)
opar<-par(mfrow=c(1,3))
svyhist(~enroll, dstrat, main="Survey weighted",col="purple",ylim=c(0,1.3e-3))
hist(apistrat$enroll,  main="Sample unweighted",col="purple",prob=TRUE,ylim=c(0,1.3e-3))
hist(apipop$enroll,  main="Population",col="purple",prob=TRUE,ylim=c(0,1.3e-3))

par(mfrow=c(1,1))
svyboxplot(enroll~stype,dstrat,all.outliers=TRUE)
svyboxplot(enroll~1,dstrat)
par(opar)
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line
