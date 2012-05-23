\name{SE}
\alias{SE}
\alias{SE.default}
\alias{SE.svrepstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract standard errors }
\description{
Extracts standard errors from an object. The default method is for
objects with a \code{\link{vcov}} method.
}
\usage{
SE(object, ...)
\method{SE}{default}(object,...)
\method{SE}{svrepstat}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object}
  \item{\dots}{Arguments for future expansion }
}
\value{
  Vector of standard errors.
}
\seealso{ \code{\link{vcov}}}

\keyword{models}% at least one, from doc/KEYWORDS

