\name{as.svrepdesign}
\alias{as.svrepdesign}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Convert a survey design to use replicate weights}
\description{
  Creates a replicate-weights survey design object from a traditional
  strata/cluster survey design object. For stratified designs these
  will be BRR weights, for unstratified designs the cluster jackkife is used.
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

}
\keyword{survey}% at least one, from doc/KEYWORDS

