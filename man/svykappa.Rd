\name{svykappa}
\alias{svykappa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cohen's kappa for agreement}
\description{
 Computes the unweighted kappa measure of agreement between two raters
 and the standard error.
}
\usage{
svykappa(formula, design, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{one-sided formula giving two measurements}
  \item{design}{survey design object}
  \item{\dots}{for future expansion }
}
\value{
  Object of class \code{svystat}
}

\seealso{ \code{\link{svycontrast}}}
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
svykappa(~comp.imp+sch.wide, dclus1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

