\name{svyby}
\alias{svyby}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Survey statistics on subsets}
\description{
Compute survey statistics on subsets of a survey defined by factors. 
}
\usage{
svyby(formula, by, design, FUN, ..., keep.var = FALSE, keep.names = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula specifying the variables to pass to \code{FUN} }
  \item{by}{A formula specifying factors that define subsets, or a list
    of factors.}
  \item{design}{A \code{svydesign} or \code{svrepdesign} object}
  \item{FUN}{A function taking a formula and survey design object as its
  first two arguments.}
  \item{\dots}{Other arguments to \code{FUN}}
  \item{keep.var}{If \code{FUN} returns a \code{svystat} object, extract
  standard errors from it}
  \item{keep.names}{Define row names based on the subsets}
}

\value{
  A data frame showing the factors and the results of \code{FUN}
}
\seealso{\code{\link{svytable}} }
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

svyby(~api99, ~stype, dclus1, svymean)
svyby(~api99, ~stype, dclus1, svyquantile, quantiles=0.5)
svyby(~api99, list(school.type=apiclus1$stype), dclus1, svymean)
svyby(~api99+api00, ~stype+sch.wide, dclus1, svymean, keep.var=TRUE)

rclus1<-as.svrepdesign(dclus1)

svyby(~api99, ~stype, rclus1, svrepmean)
svyby(~api99, ~stype, rclus1, svrepquantile, quantiles=0.5)
svyby(~api99, list(school.type=apiclus1$stype), rclus1, svrepmean)
svyby(~api99+api00, ~stype+sch.wide, rclus1, svrepmean, keep.var=TRUE)


}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
