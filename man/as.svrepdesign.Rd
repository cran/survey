\name{as.svrepdesign}
\alias{as.svrepdesign}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Convert a survey design to use replicate weights}
\description{
  Creates a replicate-weights survey design object from a traditional
  strata/cluster survey design object. 
}
\usage{
as.svrepdesign(design, type=c("auto", "JK1", "JKn", "BRR", "Fay"),
fay.rho = 0, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{design}{Object of class \code{survey.design} }
  \item{type}{Type of replicate weights. \code{"auto"} uses JKn for
    stratified, JK1 for unstratified designs}
  \item{fay.rho}{Tuning parameter for Fay's variance method }
  \item{\dots}{Other arguments to \code{\link{brrweights}}, controlling
    the handling of strata that do not have exactly 2 clusters.}
}

\value{
  Object of class \code{svrep.design}.
}

\seealso{\code{\link{brrweights}}, \code{\link{svydesign}}, \code{\link{svrepdesign}}}

\examples{
data(scd)
scddes<-svydesign(data=scd, prob=~1, id=~ambulance, strata=~ESA,
nest=TRUE, fpc=rep(5,6))
scdnofpc<-svydesign(data=scd, prob=~1, id=~ambulance, strata=~ESA,
nest=TRUE)

# convert to BRR replicate weights
scd2brr <- as.svrepdesign(scdnofpc, type="BRR")
# convert to JKn weights 
scd2jkn <- as.svrepdesign(scdnofpc, type="JKn")

# convert to JKn weights with finite population correction
scd2jknf <- as.svrepdesign(scddes, type="JKn")

svrepratio(~alive, ~arrests, design=scd2brr)
svrepratio(~alive, ~arrests, design=scd2jkn)
svrepratio(~alive, ~arrests, design=scd2jknf)

data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1<-as.svrepdesign(dclus1)

svrepmean(~api00, rclus1)
svreptotal(~enroll, rclus1)

}
\keyword{survey}% at least one, from doc/KEYWORDS

