\name{svyranktest}
\alias{svyranktest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Design-based rank tests
}
\description{
Design-based versions of two-sample rank tests.  The built-in tests are
all for location hypotheses, but the user can specify others.  
}
\usage{
svyranktest(formula, design, test = c("wilcoxon", "vanderWaerden", "median"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Model formula \code{y~g} for outcome variable \code{y} and group \code{g}
}
  \item{design}{
A survey design object
}
  \item{test}{
    Which rank test to use: Wilcoxon, van der Waerden's normal-scores
    test, Mood's test for the median, or a function \code{f(r,N)} where
    \code{r} is the rank and \code{N} the estimated population size. 
}
  \item{\dots}{ for future expansion}

}

\value{
Object of class \code{htest}
}
\references{
Lumley T, Scott AJ (2012) "Two-sample rank tests under complex
sampling". Technical Report 2012-4, Department of Statistics, University
of Auckland.

}

\seealso{
\code{\link{svyttest}}, \code{\link{svylogrank}}
}
\examples{

data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, fpc=~fpc, data=apiclus1)

svyranktest(ell~comp.imp, dclus1)
svyranktest(ell~comp.imp, dclus1, test="median")

## upper quartile
svyranktest(ell~comp.imp, dclus1, test=function(r,N) as.numeric(r>0.75*N))


quantiletest<-function(p){
	  rval<-function(r,N) as.numeric(r>(N*p))
	  attr(rval,"name")<-paste(p,"quantile")
	  rval
	}
svyranktest(ell~comp.imp, dclus1, test=quantiletest(0.5))
svyranktest(ell~comp.imp, dclus1, test=quantiletest(0.75))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{htest}% __ONLY ONE__ keyword per line
