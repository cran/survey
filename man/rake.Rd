\name{rake}
\alias{rake}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Raking of survey design}
\description{
  Raking uses iterative post-stratification to match marginal
  distributions of a survey sample to known population margins. 
}
\usage{
rake(design, sample.margins, population.margins, control = list(maxit =
10, epsilon = 1, verbose=FALSE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{A survey design with replicate weights }
  \item{sample.margins}{list of formulas or data frames describing
    sample margins}
  \item{population.margins}{list of tables describing corresponding
    population margins }
  \item{control}{\code{maxit} controls the number of
    iterations. Convergence is declared if the maximum change in a table
    entry is less than  \code{epsilon}. If \code{epsilon<1} it is
    taken to be a fraction of the total sampling weight. }
}
\details{
  
  Raking is known to converge for any table without zeros, and for any
  table with zeros for which there is a joint distribution with the
  given margins and the same pattern of zeros.   The `margins' need not
  be one-dimensional.
  
  The algorithm works by repeated calls to \code{\link{postStratify}},
  perhaps the least efficient possible implementation.
}
\value{
 A raked survey design. 
}

\seealso{\code{\link{as.svrepdesign}}, \code{\link{svrepdesign}}, \code{\link{postStratify}} }
\examples{


}
\keyword{survey}% at least one, from doc/KEYWORDS

