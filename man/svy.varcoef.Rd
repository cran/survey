\name{svy.varcoef}
\alias{svy.varcoef}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Sandwich variance estimator for glms}
\description{
Computes the sandwich variance estimator for a generalised linear model fitted to data from a complex sample survey. Designed to be used internally by \code{\link{svyglm}}.
}
\usage{
svy.varcoef(glm.object, design,std.errors = c("linearized",
                 "Bell-McCaffrey", "Bell-McCaffrey-2"), degf = FALSE)
}
\arguments{
  \item{glm.object}{A \code{\link{glm}} object}
  \item{design}{A \code{survey.design} object }
    \item{std.errors}{The kind of standard errors to compute}
  \item{degf}{Whether to compute the adjusted degrees of freedom along with 
    Bell-McCaffrey standard errors}

}
\value{
 A variance matrix
}
\author{ Thomas Lumley}


\seealso{\code{\link{svyglm}},\code{\link{svydesign}}, \code{\link{svyCprod}} }


\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
