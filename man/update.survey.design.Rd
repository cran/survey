\name{update.survey.design}
\alias{update.survey.design}
\alias{update.twophase}
\alias{update.svyrep.design}
\alias{update.DBIsvydesign}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Add variables to a survey design}
\description{
Update the data variables in a survey design, either with a formula for a new set of variables or with an expression for variables to be added.
}
\usage{
\method{update}{survey.design}(object, ...)
\method{update}{twophase}(object, ...)
\method{update}{svyrep.design}(object, ...)
\method{update}{DBIsvydesign}(object, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{a survey design object}
  \item{\dots}{Arguments \code{tag=expr} add a new variable \code{tag}
    computed by evaluating \code{expr} in the survey data.}
}
\details{
Database-backed objects may not have write access to the database and so
\code{update} does not attempt to modify the database.  The expressions
are stored and are evaluated when the data is loaded.

If a set of new variables will be used extensively it may be more efficient to
modify the database, either with SQL queries from the R interface or
separately. One useful intermediate approach is to create a table with
the new variables and a view that joins this table to the table of
existing variables.

There is now a base-R function \code{\link{transform}} for adding new
  variables to a data frame, so I have added \code{transform} as a synonym for
  \code{update} for survey objects.
  }
\value{
A survey design object
}

\seealso{\code{\link{svydesign}}, \code{\link{svrepdesign}}, \code{\link{twophase}}}

\examples{
data(api)
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat,
fpc=~fpc)
dstrat<-update(dstrat, apidiff=api00-api99)
svymean(~api99+api00+apidiff, dstrat)
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
