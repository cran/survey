\name{postStratify}
\alias{postStratify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Post-stratify a replicate weight survey }
\description{
 Post-stratification adjusts the sampling and replicate weights so that
 the joint distribution of a set of post-stratifying variables matches
 the known population joint distribution. The advantage of
 post-stratification is that sampling frames need not be available for
 the strata.
}
\usage{
postStratify(design, strata, population, partial = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{A survey design with replicate weights}
  \item{strata}{A formula or data frame of post-stratifying variables }
  \item{population}{A \code{\link{table}}, \code{\link{xtabs}} or \code{data.frame}
    with population frequencies }
  \item{partial}{if \code{TRUE}, ignore population strata not present in
  the sample}
}
\details{
  The \code{population} totals can be specified as a table with the
  strata variables in the margins, or as a data frame where one column
  lists frequencies and the other columns list the unique combinations
  of strata variables (the format produced by \code{as.data.frame}
  acting on a \code{table} object).
}
\value{
 A new survey design object.
}
\references{ }

\seealso{\code{\link{as.svrepdesign}}, \code{svrepdesign} }
\examples{

}
\keyword{survey}% at least one, from doc/KEYWORDS

