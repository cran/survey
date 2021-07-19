\name{svyttest}
\alias{svyttest}
\alias{confint.svyttest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Design-based t-test}
\description{
One-sample or two-sample t-test.  This function is a wrapper for
\code{\link{svymean}} in the one-sample case and for
\code{\link{svyglm}} in the two-sample case. Degrees of freedom are
\code{degf(design)-1} for the one-sample test and \code{degf(design)-2}
for the two-sample case.
}
\usage{
svyttest(formula, design,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Formula, \code{outcome~group} for two-sample,
    \code{outcome~0} or \code{outcome~1} for one-sample. The \code{group} variable 
    must be a factor or character with two levels, or be coded 0/1 or 1/2}
  \item{design}{survey design object}
  \item{\dots}{for methods }
}
\value{
   Object of class \code{htest}
}

\seealso{ \code{\link{t.test}}}
\examples{
data(api)
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)
tt<-svyttest(enroll~comp.imp, dclus2)
tt
confint(tt, level=0.9)

svyttest(enroll~I(stype=="E"),dclus2)

svyttest(I(api00-api99)~0, dclus2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{htest}% __ONLY ONE__ keyword per line
