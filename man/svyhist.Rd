\name{svyhist}
\alias{svyhist}
\alias{svyboxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Histograms and boxplots}
\description{
Histograms and boxplots weighted by the sampling weights. 
}
\usage{
svyhist(formula, design, xlab = NULL, main = NULL,
  probability=TRUE,freq=!probability, ...)
svyboxplot(formula,design,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{One-sided formula for \code{svyhist}, two-sided for \code{svyboxplot}}
  \item{design}{A survey design object}
  \item{xlab}{x-axis label}
  \item{main}{Main title}
  \item{probability,freq}{Y-axis is probability density or frequency}
  \item{\dots}{Other arguments to \code{\link{hist}} or \code{\link{bxp}}}
}

\details{
The histogram breakpoints are computed as if the sample were a
simple random sample of the same size.

The boxplot whiskers go to the maximum and minimum observations or to
1.5 interquartile ranges beyond the end of the box, whichever is
closer. The maximum and minimum are plotted as outliers if they are
beyond the ends of the whiskers, but other outlying points are not
plotted. This is likely to change in the future.  \code{svyboxplot}
requires a two-sided formula; use \code{variable~1} for a single boxplot.
 }

\seealso{ \code{\link{svyplot}}}
\examples{
data(api)
dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat, 
    fpc = ~fpc)
svyhist(~ell, dstrat, main="English language learners",col="peachpuff")

svyboxplot(enroll~stype,dstrat)
svyboxplot(enroll~1,dstrat)
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line
