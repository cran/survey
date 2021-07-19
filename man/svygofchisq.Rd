\name{svygofchisq}
\alias{svygofchisq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test of fit to known probabilities
}
\description{
A Rao-Scott-type version of the chi-squared test for goodness of fit to prespecified proportions.  The test statistic is the chi-squared statistic applied to the estimated population table, and the reference distribution is a Satterthwaite approximation: the test statistic divided by the estimated scale is compared to a chi-squared distribution with the estimated df.
}
\usage{
svygofchisq(formula, p, design, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Formula specifying a single factor variable
}
  \item{p}{
Vector of probabilities for the categories of the factor, in the correct order (will be rescaled to sum to 1)
}
  \item{design}{
Survey design object
}
  \item{\dots}{
Other arguments to pass to \code{\link{svytotal}}, such as \code{na.rm}
}
}

\value{
An object of class \code{htest}
}

\seealso{
\code{\link{chisq.test}}, \code{\link{svychisq}}, \code{\link{pchisqsum}}
}
\examples{
data(api)
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)

true_p <- table(apipop$stype)

svygofchisq(~stype,dclus2,p=true_p)
svygofchisq(~stype,dclus2,p=c(1/3,1/3,1/3))

}
