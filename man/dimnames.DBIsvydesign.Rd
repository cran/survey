\name{dimnames.DBIsvydesign}
\alias{dimnames.DBIsvydesign}
\alias{dimnames.ODBCsvydesign}
\alias{dimnames.survey.design}
\alias{dimnames.svyrep.design}
\alias{dimnames.twophase}
\alias{dimnames.svyimputationList}
\alias{dim.DBIsvydesign}
\alias{dim.ODBCsvydesign}
\alias{dim.survey.design}
\alias{dim.twophase}
\alias{dim.svyimputationList}
\alias{dim.svyrep.design}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dimensions of survey designs}
\description{
\code{dimnames} returns variable names and row names for the data
variables in a design object and \code{dim} returns dimensions.
For multiple imputation designs there is a third dimension giving the
number of imputations.  For database-backed designs the second dimension
includes variables defined by \code{update}.  The first dimension
excludes observations with zero weight.
}
\usage{
\method{dim}{survey.design}(x)
\method{dim}{svyimputationList}(x)
\method{dimnames}{survey.design}(x)
\method{dimnames}{DBIsvydesign}(x)
\method{dimnames}{ODBCsvydesign}(x)
\method{dimnames}{svyimputationList}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design object}

}

\value{
A vector of numbers for \code{dim}, a list of vectors of strings for \code{dimnames}.
}

\seealso{ \code{\link{update.DBIsvydesign}}, \code{\link{with.svyimputationList}}}
\examples{
data(api)
dclus1 <- svydesign(ids=~dnum,weights=~pw,data=apiclus1,fpc=~fpc)
dim(dclus1)
dimnames(dclus1)
colnames(dclus1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{manip}% __ONLY ONE__ keyword per line
