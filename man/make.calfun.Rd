\name{make.calfun}
\alias{make.calfun}
\alias{cal.linear}
\alias{cal.raking}
\alias{cal.logit}
\alias{cal.sinh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calibration metrics}
\description{
Create calibration metric for use  in \code{\link{calibrate}}. The
function \code{F} is the link function described in section 2 of
Deville et al. To create a new calibration metric, specify \eqn{F-1}{F-1} and its
derivative. The package provides \code{cal.linear}, \code{cal.raking},
\code{cal.logit}, which are standard, and \code{cal.sinh} from the
\code{CALMAR2} macro, for which \code{F} is the derivative of the inverse hyperbolic
sine. 
}
\usage{
make.calfun(Fm1, dF, name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fm1}{Function \eqn{F-1}{F-1} taking a vector \code{u} and a
    vector of length 2, \code{bounds}.}
  \item{dF}{Derivative of \code{Fm1} wrt \code{u}: arguments \code{u}
    and \code{bounds} }
  \item{name}{Character string to use as name }
}
\value{
An object of class \code{"calfun"}
}

\references{
Deville J-C, Sarndal C-E, Sautory O (1993) Generalized Raking
Procedures in Survey Sampling. JASA 88:1013-1020

Deville J-C, Sarndal C-E (1992) Calibration Estimators in Survey
Sampling. JASA 87: 376-382
 }

\seealso{\code{\link{calibrate}} }
\examples{
str(cal.linear)
cal.linear$Fm1
cal.linear$dF

hellinger <- make.calfun(Fm1=function(u, bounds)  ((1-u/2)^-2)-1,
                    dF= function(u, bounds) (1-u/2)^-3 ,
                    name="hellinger distance")

hellinger

data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

svymean(~api00,calibrate(dclus1, ~api99, pop=c(6194, 3914069),
         calfun=hellinger))

svymean(~api00,calibrate(dclus1, ~api99, pop=c(6194, 3914069),
         calfun=cal.linear))

svymean(~api00,calibrate(dclus1, ~api99, pop=c(6194,3914069),
          calfun=cal.raking))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

