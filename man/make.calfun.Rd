\name{make.calfun}
\alias{cal.linear}
\alias{cal.raking}
\alias{cal.logit}
\alias{make.calfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calibration metrics}
\description{
\code{make.calfun()} creates objects representing distance functions for
calibration of weights.  Three such functions are pre-defined.
}
\usage{
make.calfun(Fm1, dF, name)
cal.linear(u, bounds)
cal.raking(u, bounds)
cal.logit(u, bounds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fm1}{Link function minus 1 for calibration factors, function of two
    arguments: \code{u} and \code{bounds} }
  \item{dF}{Derivative of \code{Fm1} with respect to \code{u} }
  \item{name}{Character string giving a name for the object}
  \item{u}{Linear predictor}
  \item{bounds}{\code{c(lower,upper)} bounds for calibration factors}
}

\value{
For \code{make.calfun}, an object of class \code{calfun}. 
}
\references{
  Deville J-C, Sarndal C-E. (1992) Calibration Estimators in Survey
  Sampling. JASA 87:376-382
  }
\note{
The link function F is the same as F in Table 1 of Deville and Sarndal (1992)
}
\seealso{\code{\link{calibrate}} }
\examples{

hellinger<-make.calfun(function(u, bounds) ((1-u/2)^-2)-1,
                       function(u, bounds) (1-u/2)^{-3},
                       "hellinger distance")
hellinger

data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)


svymean(~api00,calibrate(dclus1, ~api99, pop=c(6194, 3914069),calfun=hellinger))

svymean(~api00,calibrate(dclus1, ~api99, pop=c(6194, 3914069),calfun=cal.linear))

svymean(~api00,calibrate(dclus1, ~api99, pop=c(6194, 3914069),calfun=cal.raking))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

