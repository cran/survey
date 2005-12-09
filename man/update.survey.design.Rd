\name{update.survey.design}
\alias{update.survey.design}
\alias{update.twophase}
\alias{update.svyrep.design}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Add variables to a survey design}
\description{
Update the data variables in a survey design, either with a formula for a new set of variables or with an expression for variables to be added.
}
\usage{
\method{update}{survey.design}(object, ...)
\method{update}{twophase}(object, ...)
\method{update}{svyrep.design}(object, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{a survey design object}
  \item{\dots}{Arguments \code{tag=expr} add a new variable \code{tag}
    computed by evaluating \code{expr} in the survey data.}
}

\value{
A survey design object
}

\seealso{\code{\link{svydesign}}, \code{\link{svrepdesign}}, \code{\link{twophase}}}

\examples{
data(api)
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat,
fpc=~fpc)
dstrat<-update(dstrat, apidiff=api00-api99)
svymean(~api99+api00+apidiff, dstrat)
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
