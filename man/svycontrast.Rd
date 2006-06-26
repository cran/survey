\name{svycontrast}
\alias{svycontrast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Linear constrasts of survey statistics }
\description{
Computes linear contrasts from objects of class \code{svrepstat} or \code{svystat}.
}
\usage{
svycontrast(stat, contrasts, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stat}{object of class \code{svrepstat} or \code{svystat} }
  \item{contrasts}{A vector or list of vectors of contrasts }
  \item{\dots}{For future expansion}
}
\value{
Object of class \code{svrepstat} or \code{svystat}
}

\seealso{\code{\link{regTermTest}}, \code{svyglm}}
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

a <- svytotal(~api00+api99, dclus1)
svycontrast(a, list(avg=c(0.5,0.5), diff=c(1,-1)))

}
\keyword{survey}% at least one, from doc/KEYWORDS

