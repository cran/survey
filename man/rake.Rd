\name{rake}
\alias{rake}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Raking of replicate weight design}
\description{
  Raking uses iterative post-stratification to match marginal
  distributions of a survey sample to known population margins. 
}
\usage{
rake(design, sample.margins, population.margins, control = list(maxit =
10, epsilon = 1, verbose=FALSE), compress=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{A survey object }
  \item{sample.margins}{list of formulas or data frames describing
    sample margins}
  \item{population.margins}{list of tables or data frames
    describing corresponding population margins }
  \item{control}{\code{maxit} controls the number of
    iterations. Convergence is declared if the maximum change in a table
    entry is less than  \code{epsilon}. If \code{epsilon<1} it is
    taken to be a fraction of the total sampling weight. }
  \item{compress}{If \code{design} has replicate weights, attempt to
    compress the new replicate weight matrix? When \code{NULL}, will
    attempt to compress if the original weight matrix was compressed}
}
\details{
  The \code{sample.margins} should be in a format suitable for \code{\link{postStratify}}.
  
  Raking (aka iterative proportional fitting) is known to converge for
  any table without zeros, and for any table with zeros for which there
  is a joint distribution with the given margins and the same pattern of
  zeros.  The `margins' need not be one-dimensional.
  
  The algorithm works by repeated calls to \code{\link{postStratify}},
  not the most efficient possible implementation.
}
\value{
 A raked survey design. 
}

\seealso{\code{\link{as.svrepdesign}}, \code{\link{svrepdesign}},
  \code{\link{postStratify}}, \code{\link{compressWeights}} }
\examples{
data(api)
dclus1 <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1 <- as.svrepdesign(dclus1)

svrepmean(~api00, rclus1)
svreptotal(~enroll, rclus1)

## population marginal totals for each stratum
pop.types <- data.frame(stype=c("E","H","M"), Freq=c(4421,755,1018))
pop.schwide <- data.frame(sch.wide=c("No","Yes"), Freq=c(1072,5122))

rclus1r <- rake(rclus1, list(~stype,~sch.wide), list(pop.types, pop.schwide))

svrepmean(~api00, rclus1r)
svreptotal(~enroll, rclus1r)

## marginal totals correspond to population
xtabs(~stype, apipop)
svreptable(~stype, rclus1r, round=TRUE)
xtabs(~sch.wide, apipop)
svreptable(~sch.wide, rclus1r, round=TRUE)

## joint totals don't correspond 
xtabs(~stype+sch.wide, apipop)
svreptable(~stype+sch.wide, rclus1r, round=TRUE)

## compare to joint post-stratification
## (only possible if joint population table is known)
##
pop.table <- xtabs(~stype+sch.wide,apipop)
rclus1ps <- postStratify(rclus1, ~stype+sch.wide, pop.table)
svreptable(~stype+sch.wide, rclus1ps, round=TRUE)

svrepmean(~api00, rclus1ps)
svreptotal(~enroll, rclus1ps)

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}
