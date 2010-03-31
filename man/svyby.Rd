\name{svyby}
\alias{svyby}
\alias{svyby.default}
\alias{SE.svyby}
\alias{deff.svyby}
\alias{coef.svyby}
\alias{unwtd.count}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Survey statistics on subsets}
\description{
Compute survey statistics on subsets of a survey defined by factors. 
}
\usage{
svyby(formula, by ,design,...)
\method{svyby}{default}(formula, by, design, FUN, ..., deff=FALSE,keep.var = TRUE,
keep.names = TRUE,verbose=FALSE, vartype=c("se","ci","ci","cv","cvpct","var"),
 drop.empty.groups=TRUE, covmat=FALSE, return.replicates=FALSE, multicore=getOption("survey.multicore"))
\method{SE}{svyby}(object,...)
\method{deff}{svyby}(object,...)
\method{coef}{svyby}(object,...)
unwtd.count(x, design, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula,x}{A formula specifying the variables to pass to
    \code{FUN} (or a matrix, data frame, or vector)}
  \item{by}{A formula specifying factors that define subsets, or a list
    of factors.}
  \item{design}{A \code{svydesign} or \code{svrepdesign} object}
  \item{FUN}{A function taking a formula and survey design object as its
  first two arguments.}
  \item{\dots}{Other arguments to \code{FUN}}
  \item{deff}{Request a design effect from \code{FUN}}
  \item{keep.var}{If \code{FUN} returns a \code{svystat} object, extract
  standard errors from it}
\item{keep.names}{Define row names based on the subsets}
\item{verbose}{If \code{TRUE}, print a label for each subset as it is
  processed.}
\item{vartype}{Report variability as one or more of
  standard error, confidence interval, coefficient of
  variation,  percent coefficient of variation, or variance}
\item{drop.empty.groups}{If \code{FALSE}, report \code{NA} for empty
  groups, if \code{TRUE} drop them from the output}
\item{covmat}{If \code{TRUE}, compute covariances between estimates for
  different subsets (currently only for replicate-weight
  designs). Allows \code{\link{svycontrast}} to be used on output.}
\item{return.replicates}{Only for replicate-weight designs. If
  \code{TRUE}, return all the replicates as the "replicates" attribute of the result}
\item{multicore}{Use \code{multicore} package to distribute subsets over
   multiple processors?}
\item{object}{An object of class \code{"svyby"}}
}

\value{
  An object of class \code{"svyby"}: a data frame showing the factors and the results of \code{FUN}.

For \code{unwtd.count}, the unweighted number of non-missing observations in the data matrix specified by \code{x} for the design. 
}
\details{
The variance type "ci" asks for confidence intervals, which are produced
  by \code{confint}. In some cases additional options to \code{FUN} will
  be needed to produce confidence intervals, for example,
  \code{svyquantile} needs \code{ci=TRUE} or \code{keep.var=FALSE}.

  \code{unwtd.count} is designed to be passed to \code{svyby} to report
  the number of non-missing observations in each subset. Observations
  with exactly zero weight will also be counted as missing, since that's
  how subsets are implemented for some designs.

  Parallel processing with \code{multicore=TRUE} is useful only for
  fairly large problems and on computers with sufficient memory. The
  \code{multicore} package is incompatible with some GUIs, although the
  Mac Aqua GUI appears to be safe. 

  }

\note{ Asking for a design effect (\code{deff=TRUE}) from a function
  that does not produce one will cause an error or incorrect formatting
  of the output. The same will occur with \code{keep.var=TRUE} if the
  function does not compute a standard error.
}
\seealso{\code{\link{svytable}} and \code{\link{ftable.svystat}} for
  contingency tables, \code{\link{ftable.svyby}} for pretty-printing of \code{svyby} }
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

svyby(~api99, ~stype, dclus1, svymean)
svyby(~api99, ~stype, dclus1, svyquantile, quantiles=0.5,ci=TRUE,vartype="ci")
## without ci=TRUE svyquantile does not compute standard errors
svyby(~api99, ~stype, dclus1, svyquantile, quantiles=0.5, keep.var=FALSE)
svyby(~api99, list(school.type=apiclus1$stype), dclus1, svymean)
svyby(~api99+api00, ~stype, dclus1, svymean, deff=TRUE,vartype="ci")
svyby(~api99+api00, ~stype+sch.wide, dclus1, svymean, keep.var=FALSE)
## report raw number of observations
svyby(~api99+api00, ~stype+sch.wide, dclus1, unwtd.count, keep.var=FALSE)

rclus1<-as.svrepdesign(dclus1)

svyby(~api99, ~stype, rclus1, svymean)
svyby(~api99, ~stype, rclus1, svyquantile, quantiles=0.5)
svyby(~api99, list(school.type=apiclus1$stype), rclus1, svymean, vartype="cv")
svyby(~enroll,~stype, rclus1,svytotal, deff=TRUE)
svyby(~api99+api00, ~stype+sch.wide, rclus1, svymean, keep.var=FALSE)
##report raw number of observations
svyby(~api99+api00, ~stype+sch.wide, rclus1, unwtd.count, keep.var=FALSE)

## comparing subgroups using covmat=TRUE
mns<-svyby(~api99, ~stype, rclus1, svymean,covmat=TRUE)
vcov(mns)
svycontrast(mns, c(E = 1, M = -1))

str(svyby(~api99, ~stype, rclus1, svymean,return.replicates=TRUE))


## extractor functions
(a<-svyby(~enroll, ~stype, rclus1, svytotal, deff=TRUE, verbose=TRUE, vartype=c("se","cv","cvpct","var")))
deff(a)
SE(a)
cv(a)
coef(a)

## ratio estimates
svyby(~api.stu, by=~stype, denominator=~enroll, design=dclus1, svyratio)

## empty groups
svyby(~api00,~comp.imp+sch.wide,design=dclus1,svymean)
svyby(~api00,~comp.imp+sch.wide,design=dclus1,svymean,drop.empty.groups=FALSE)

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
