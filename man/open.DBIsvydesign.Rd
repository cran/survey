\name{open.DBIsvydesign}
\alias{open.DBIsvydesign}
\alias{close.DBIsvydesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Open and close DBI connections }
\description{
  A database-backed survey design object contains a connection to a
database.  This connection will be broken if the object is saved and
reloaded, and the connection should ideally be closed with \code{close}
before quitting R (although it doesn't matter for SQLite
connections). The connection can be reopened with \code{open}.
}
\usage{
\method{open}{DBIsvydesign}(con, ...)
\method{close}{DBIsvydesign}(con, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{Object of class \code{DBIsvydesign} }
  \item{\dots}{Other options, to be passed to \code{dbConnect} or
    \code{dbDisconnect}.}
}
\value{
The same survey design object with the connection opened or closed.
}

\seealso{\code{\link{svydesign}}

  DBI package }
\examples{
\dontrun{
library(RSQLite)
dbclus1<-svydesign(id=~dnum, weights=~pw, fpc=~fpc,
data="apiclus1",dbtype="SQLite",
dbname=system.file("api.db",package="survey"))

dbclus1
close(dbclus1)
dbclus1
try(svymean(~api00, dbclus1))

dbclus1<-open(dbclus1)
open(dbclus1)
svymean(~api00, dbclus1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

