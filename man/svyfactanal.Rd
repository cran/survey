\name{svyfactanal}
\Rdversion{1.1}
\alias{svyfactanal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Factor analysis in complex surveys (experimental).
}
\description{
This function fits a factor analysis model or SEM,  by maximum weighted likelihood.
}
\usage{
svyfactanal(formula, design, factors, 
   n = c("none", "sample", "degf","effective", "min.effective"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Model formula specifying the variables to use
}
  \item{design}{
Survey design object
}
  \item{factors}{
Number of factors to estimate
}
  \item{n}{
Sample size to be used for testing: see below}
  \item{\dots}{
Other arguments to pass to \code{\link{factanal}}.
}
}

\details{
The population covariance matrix is estimated by \code{\link{svyvar}}
and passed to \code{\link{factanal}}

Although fitting these models requires only the estimated covariance
matrix, inference requires a sample size.  With \code{n="sample"}, the sample size is taken to be
the number of observations; with \code{n="degf"}, the survey degrees of
freedom as returned by \code{\link{degf}}.  Using \code{"sample"}
corresponds to standardizing weights to have mean 1, and is known to
result in anti-conservative tests.

The other two methods estimate an effective sample size for each
variable as the sample size where the standard error of a variance of a
Normal distribution would match the design-based standard error
estimated by \code{\link{svyvar}}. With \code{n="min.effective"} the
minimum sample size across the variables is used; with
\code{n="effective"} the harmonic mean is used.  For \code{svyfactanal}
the test of model adequacy is optional, and the default choice,
\code{n="none"}, does not do the test.

}
\value{
An object of class \code{factanal}
}
\references{
.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{factanal}}

The \code{lavaan.survey} package fits structural equation models to complex samples using similar techniques.
}
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

svyfactanal(~api99+api00+hsg+meals+ell+emer, design=dclus1, factors=2)

svyfactanal(~api99+api00+hsg+meals+ell+emer, design=dclus1, factors=2, n="effective")

##Population dat for comparison
factanal(~api99+api00+hsg+meals+ell+emer, data=apipop, factors=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{multivariate}% __ONLY ONE__ keyword per line
