\name{svrVar}
\alias{svrVar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Compute variance from replicates }
\description{
Compute an appropriately scaled empirical variance estimate from
replicates.  The \code{mse} argument specifies whether the sums of
squares should be centered at the point estimate (\code{mse=TRUE}) or
the mean of the replicates. It is usually taken from the \code{mse}
component of the design object.
}
\usage{
svrVar(thetas, scale, rscales, na.action=getOption("na.action"), mse=getOption("survey.replicates.mse"),coef)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{thetas}{matrix whose rows are replicates (or a vector of replicates)}
  \item{scale}{Overall scaling factor}
  \item{rscales}{Scaling factor for each squared deviation }
  \item{na.action}{How to handle replicates where the statistic could
    not be estimated}
  \item{mse}{if \code{TRUE}, center at the point estimated, if
    \code{FALSE} center at the mean of the replicates}
  \item{coef}{The point estimate, required only if \code{mse==TRUE}}
}
\value{
covariance matrix.  
}
\seealso{\code{\link{svrepdesign}}, \code{\link{as.svrepdesign}},
  \code{\link{brrweights}},
  \code{\link{jk1weights}}, \code{\link{jknweights}}}

\keyword{survey}% at least one, from doc/KEYWORDS

