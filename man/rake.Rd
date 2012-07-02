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
    sample margins, which must not contain missing values}
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
  
  The algorithm works by repeated calls to \code{\link{postStratify}}
  (iterative proportional fitting), which is efficient for large
  multiway tables. For small tables \code{\link{calibrate}} will be
  faster, and also allows raking to population totals for continuous
  variables, and raking with bounded weights. 
}
\value{
 A raked survey design. 
}

\seealso{
  \code{\link{postStratify}}, \code{\link{compressWeights}}
  
  \code{\link{calibrate}} for other ways to use auxiliary information.
}
\examples{
data(api)
dclus1 <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1 <- as.svrepdesign(dclus1)

svymean(~api00, rclus1)
svytotal(~enroll, rclus1)

## population marginal totals for each stratum
pop.types <- data.frame(stype=c("E","H","M"), Freq=c(4421,755,1018))
pop.schwide <- data.frame(sch.wide=c("No","Yes"), Freq=c(1072,5122))

rclus1r <- rake(rclus1, list(~stype,~sch.wide), list(pop.types, pop.schwide))

svymean(~api00, rclus1r)
svytotal(~enroll, rclus1r)

## marginal totals correspond to population
xtabs(~stype, apipop)
svytable(~stype, rclus1r, round=TRUE)
xtabs(~sch.wide, apipop)
svytable(~sch.wide, rclus1r, round=TRUE)

## joint totals don't correspond 
xtabs(~stype+sch.wide, apipop)
svytable(~stype+sch.wide, rclus1r, round=TRUE)

## Do it for a design without replicate weights
dclus1r<-rake(dclus1, list(~stype,~sch.wide), list(pop.types, pop.schwide))

svymean(~api00, dclus1r)
svytotal(~enroll, dclus1r)

## compare to raking with calibrate()
dclus1gr<-calibrate(dclus1, ~stype+sch.wide, pop=c(6194, 755,1018,5122),
           calfun="raking")
svymean(~stype+api00, dclus1r)
svymean(~stype+api00, dclus1gr)

## compare to joint post-stratification
## (only possible if joint population table is known)
##
pop.table <- xtabs(~stype+sch.wide,apipop)
rclus1ps <- postStratify(rclus1, ~stype+sch.wide, pop.table)
svytable(~stype+sch.wide, rclus1ps, round=TRUE)

svymean(~api00, rclus1ps)
svytotal(~enroll, rclus1ps)

## Example of raking with partial joint distributions
pop.imp<-data.frame(comp.imp=c("No","Yes"),Freq=c(1712,4482))
dclus1r2<-rake(dclus1, list(~stype+sch.wide, ~comp.imp),
               list(pop.table, pop.imp))
svymean(~api00, dclus1r2)

## compare to calibrate() syntax with tables
dclus1r2<-calibrate(dclus1, formula=list(~stype+sch.wide, ~comp.imp),
               population=list(pop.table, pop.imp),calfun="raking")
svymean(~api00, dclus1r2)


}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}
