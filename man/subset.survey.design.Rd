\name{subset.survey.design}
\alias{subset.survey.design}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Subset of survey}
\description{
Restrict a survey design to a subpopulation, keeping the original design information about number of clusters, strata. 
}
\usage{
\method{subset}{survey.design}(x, subset, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A \code{survey.design} object}
  \item{subset}{An expression specifying the subpopulation}
  \item{\dots}{Arguments not used by this method}
}
\value{
  A \code{survey.design} object
}

\seealso{\code{\link{svydesign}}}

\examples{
data(fpc)
dfpc<-svydesign(id=~psuid,strat=~stratid,weight=~weight,data=fpc,nest=TRUE)
dsub<-subset(dfpc,x>4)
summary(dsub)
## These should give the same estimates and variances
svymean(~x,dsub)
svyglm(x~I(x>4)+0,design=dfpc)
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
