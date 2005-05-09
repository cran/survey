\name{subset.survey.design}
\alias{subset.survey.design}
\alias{subset.svyrep.design}
\alias{[.survey.design}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Subset of survey}
\description{
Restrict a survey design to a subpopulation, keeping the original design
information about number of clusters, strata.   If the design has no
post-stratification or calibration data the subset will use
proportionately less memory.
}
\usage{
\method{subset}{survey.design}(x, subset, ...)
\method{subset}{svyrep.design}(x, subset, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A survey design object}
  \item{subset}{An expression specifying the subpopulation}
  \item{\dots}{Arguments not used by this method}
}
\value{
  A new survey design object
}

\seealso{\code{\link{svydesign}}}

\examples{
data(fpc)
dfpc<-svydesign(id=~psuid,strat=~stratid,weight=~weight,data=fpc,nest=TRUE)
dsub<-subset(dfpc,x>4)
summary(dsub)
svymean(~x,design=dsub)

## These should give the same domain estimates and standard errors
svyby(~x,~I(x>4),design=dfpc, svymean)
summary(svyglm(x~I(x>4)+0,design=dfpc))

data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1<-as.svrepdesign(dclus1)
svymean(~enroll, subset(dclus1, sch.wide=="Yes" & comp.imp=="Yes"))
svymean(~enroll, subset(rclus1, sch.wide=="Yes" & comp.imp=="Yes"))

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
