\name{svycdf}
\alias{svycdf}
\alias{print.svycdf}
\alias{plot.svycdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cumulative Distribution Function}
\description{
Estimates the population cumulative distribution function for specified
variables.  In contrast to \code{\link{svyquantile}}, this does not do
any interpolation: the result is a right-continuous step function.
}
\usage{
svycdf(formula, design, na.rm = TRUE,...)
\method{print}{svycdf}(x,...)
\method{plot}{svycdf}(x,xlab=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{one-sided formula giving variables from the design object }
  \item{design}{survey design object }
  \item{na.rm}{remove missing data (case-wise deletion)?}
  \item{...}{other arguments to \code{\link{plot.stepfun}}}
  \item{x}{object of class \code{svycdf}}
  \item{xlab}{a vector of x-axis labels or \code{NULL} for the default labels}
}

\value{
 An object of class \code{svycdf}, which is a list of step functions (of
 class \code{\link{stepfun}})
 }


\seealso{ \code{\link{svyquantile}}, \code{\link{svyhist}}, \code{\link{plot.stepfun}}}
\examples{
data(api)
dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat, 
    fpc = ~fpc)
cdf.est<-svycdf(~enroll+api00+api99, dstrat)
cdf.est
## function
cdf.est[[1]]
## evaluate the function
cdf.est[[1]](800)
cdf.est[[2]](800)

## compare to population and sample CDFs.
opar<-par(mfrow=c(2,1))
cdf.pop<-ecdf(apipop$enroll)
cdf.samp<-ecdf(apistrat$enroll)
plot(cdf.pop,main="Population vs sample", xlab="Enrollment")
lines(cdf.samp,col.points="red")

plot(cdf.pop, main="Population vs estimate", xlab="Enrollment")
lines(cdf.est[[1]],col.points="red")

par(opar)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{hplot}% __ONLY ONE__ keyword per line
