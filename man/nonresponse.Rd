\name{nonresponse}
\alias{nonresponse}
\alias{sparseCells}
\alias{neighbours}
\alias{joinCells}
\alias{weights.nonresponse}
\alias{print.nonresponse}
\alias{print.nonresponseSubset}
\alias{[.nonresponse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Experimental: Construct non-response weights}
\description{
  Functions to simplify the construction of non-reponse weights by
  combining strata with small numbers or large weights.
}
\usage{
nonresponse(sample.weights, sample.counts, population)
sparseCells(object, count=0,totalweight=Inf, nrweight=1.5)
neighbours(index,object)
joinCells(object,a,...)
\method{weights}{nonresponse}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.weights}{table of sampling weight by stratifying variables}
  \item{sample.counts}{table of sample counts by stratifying variables}
  \item{population}{table of population size by stratifying variables}
  \item{object}{object of class \code{"nonresponse"}}
  \item{count}{Cells with fewer sampled units than this are "sparse"}
  \item{nrweight}{Cells with higher non-response weight than this are "sparse"}
  \item{totalweight}{Cells with average sampling weight times
    non-response weight higher than this are "sparse"}
  \item{index}{Number of a cell whose neighbours are to be found}
  \item{a,...}{Cells to join}
}
\details{
  
  When a stratified survey is conducted with imperfect response it is
  desirable to rescale the sampling weights to reflect the nonresponse.
  If some strata have small sample size, high non-response, or already
  had high sampling weights it may be desirable to get less variable
  non-response weights by averaging non-response across strata.
  Suitable strata to collapse may be similar on the stratifying
  variables and/or on the level of non-response.
  
  \code{nonresponse()} combines stratified tables of population size,
  sample size, and sample weight into an object. \code{sparseCells}
  identifies cells that may need combining. \code{neighbours} describes the
  cells adjacent to a specified cell, and \code{joinCells} collapses
  the specified cells.  When the collapsing is complete, use
  \code{weights()} to extract the nonresponse weights.
  
  
}
\value{
  \code{nonresponse} and \code{joinCells} return objects of class \code{"nonresponse"},
  \code{neighbours} and \code{sparseCells} return objects of class \code{"nonresponseSubset"}
}

\examples{
data(api)
## pretend the sampling was stratified on three variables
poptable<-xtabs(~sch.wide+comp.imp+stype,data=apipop)
sample.count<-xtabs(~sch.wide+comp.imp+stype,data=apiclus1)
sample.weight<-xtabs(pw~sch.wide+comp.imp+stype, data=apiclus1)

## create a nonresponse object
nr<-nonresponse(sample.weight,sample.count, poptable)

## sparse cells
sparseCells(nr)

## Look at neighbours
neighbours(3,nr)
neighbours(11,nr)

## Collapse some contiguous cells
nr1<-joinCells(nr,3,5,7)

## sparse cells now
sparseCells(nr1)
nr2<-joinCells(nr1,3,11,8)

nr2

## one relatively sparse cell
sparseCells(nr2)
## but nothing suitable to join it to
neighbours(3,nr2)

## extract the weights
weights(nr2)
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
