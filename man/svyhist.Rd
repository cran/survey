\name{svyhist}
\alias{svyhist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Histograms}
\description{
Draws a histogram weighted by the sampling weights.
}
\usage{
svyhist(formula, design, xlab = NULL, main = NULL,
  probability=TRUE,freq=!probability, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{One-sided formula giving the variable to plot}
  \item{design}{A survey design object}
  \item{xlab}{x-axis label}
  \item{main}{Main title}
  \item{probability,freq}{Y-axis is probability density or frequency}
  \item{\dots}{Other arguments to \code{\link{hist}}}
}

\seealso{ \code{\link{svyplot}}}
\examples{
data(api)
dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat, 
    fpc = ~fpc)
svyhist(~ell, dstrat, main="English language learners",col="peachpuff")
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line
