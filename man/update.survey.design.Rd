\name{update.survey.design}
\alias{update.survey.design}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Add variables to a survey design}
\description{
Update the data variables in a survey design, either with a formula for a new set of variables or with an expression for variables to be added.
}
\usage{
update.survey.design(object, vars = ~., ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{a \code{survey.design} object}
  \item{vars}{A formula or quoted expression. See Details  below}
  \item{\dots}{Other variables not used by this method}
}
\details{
Variables in the formula or expression are looked up first in the current variables data frame and then in the calling frame.

A formula should typically have a \code{.} on the right-hand side to indicate that the current variables are to be retained.  When \code{vars} is an expression the current variables are automatically retained.
}
\value{
A \code{survey.design} object  
}

\seealso{\code{\link{svydesign}}}

\examples{

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
