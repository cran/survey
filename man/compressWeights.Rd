\name{compressWeights}
\alias{compressWeights}
\alias{compressWeights.default}
\alias{compressWeights.repweights_compressed}
\alias{[.repweights_compressed}
\alias{dim.repweights_compressed}
\alias{dimnames.repweights_compressed}
\alias{as.matrix.repweights_compressed}
\alias{as.matrix.repweights}
\alias{as.vector.repweights_compressed}
\alias{compressWeights.svyrep.design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compress replicate weight matrix}
\description{
  Many replicate weight matrices have redundant rows, such as when
  weights are the same for all observations in a PSU.  This function
  produces a compressed form. Methods for \code{as.matrix} and
  \code{as.vector} extract and expand the weights.
}
\usage{
compressWeights(rw, ...)
\method{compressWeights}{svyrep.design}(rw,...)
\method{as.matrix}{repweights_compressed}(x,...)
\method{as.vector}{repweights_compressed}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rw}{A set of replicate weights or a \code{svyrep.design} object}
  \item{x}{A compressed set of replicate weights}
  \item{\dots}{For future expansion}
}

\value{
  An object of class \code{repweights_compressed} or a
  \code{svyrep.design} object with \code{repweights} element of class \code{repweights_compressed}
}
\seealso{\code{\link{jknweights}},\code{\link{as.svrepdesign}}}
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1c<-as.svrepdesign(dclus1,compress=TRUE)
rclus1<-as.svrepdesign(dclus1,compress=FALSE)
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
