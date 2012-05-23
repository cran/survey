\name{as.svydesign2}
\alias{as.svydesign2}
\alias{.svycheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Update to the new survey design format}
\description{
The structure of survey design objects changed in version 2.9, to allow
standard errors based on multistage sampling. \code{as.svydesign} converts an
object to the new structure and  \code{.svycheck} warns if an object
does not have the new structure.

You can set \code{options(survey.want.obsolete=TRUE)} to suppress the
warnings produced by \code{.svycheck} and
\code{options(survey.ultimate.cluster=TRUE)} to always compute
variances based on just the first stage of sampling.
}
\usage{
as.svydesign2(object)
.svycheck(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{produced by \code{svydesign}}
}

\value{
Object of class \code{survey.design2}
}

\seealso{\code{\link{svydesign}}, \code{\link{svyrecvar}}}

\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
