\name{svy.varcoef}
\alias{svy.varcoef}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Sandwich variance estimator for glms}
\description{
Computes the sandwich variance estimator for a generalised linear model fitted to data from a complex sample survey. Designed to be used internally by \code{\link{svyglm}}.
}
\usage{
svy.varcoef(glm.object, design)
}
\arguments{
  \item{glm.object}{A \code{\link{glm}} object}
  \item{design}{A \code{survey.design} object }
}
\value{
 A variance matrix
}
\author{ Thomas Lumley}


\seealso{\code{\link{svyglm}},\code{\link{svydesign}}, \code{\link{svyCprod}} }


\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
