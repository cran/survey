\name{weights.survey.design}
\alias{weights.survey.design}
\alias{weights.svyrep.design}
\alias{weights.survey_fpc}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Survey design weights}
\description{
Extract weights from a survey design object.
}
\usage{
\method{weights}{survey.design}(object, ...)
\method{weights}{svyrep.design}(object,
type=c("replication","sampling","analysis"), ...)
\method{weights}{survey_fpc}(object,final=TRUE,...)
} 
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{Survey design object}
  \item{type}{Type of weights: \code{"analysis"} combines sampling and
    replication weights.}
  \item{final}{If \code{FALSE} return a data frame with sampling
    weights at each stage of sampling.}
  \item{\dots}{Other arguments ignored }
}

\value{
  vector or matrix of weights
}

\seealso{\code{\link{svydesign}}, \code{\link{svrepdesign}},
  \code{\link{as.fpc}} }

\examples{
data(scd)


scddes<-svydesign(data=scd, prob=~1, id=~ambulance, strata=~ESA,
                 nest=TRUE, fpc=rep(5,6))
repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1), c(0,1,0,1,1,0))
scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)

weights(scdrep)
weights(scdrep, type="sampling")
weights(scdrep, type="analysis")
weights(scddes)

}
\keyword{survey}% at least one, from doc/KEYWORDS

