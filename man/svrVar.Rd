\name{svrVar}
\alias{svrVar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Compute variance from replicates }
\description{
Compute an appropriately scaled empirical variance estimate from
replicates.  
}
\usage{
svrVar(thetas, scale, rscales, na.action=getOption("na.action"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{thetas}{matrix whose rows are replicates (or a vector of replicates)}
  \item{scale}{Overall scaling factor}
  \item{rscales}{Scaling factor for each squared deviation }
  \item{na.action}{How to handle replicates where the statistic could
    not be estimated}
}
\value{
covariance matrix.  
}
\seealso{\code{\link{svrepdesign}}, \code{\link{as.svrepdesign}},
  \code{\link{brrweights}},
  \code{\link{jk1weights}}, \code{\link{jknweights}}}
\examples{

}
\keyword{survey}% at least one, from doc/KEYWORDS

