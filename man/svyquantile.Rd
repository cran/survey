\name{svyquantile}
\alias{svyquantile}
\alias{svytable}
\alias{svreptable}
\alias{svymean}
\alias{svrepmean}
\alias{svytotal}
\alias{print.svystat}
\alias{print.svrepstat}
\alias{svreptotal}
\alias{svyvar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summary statistics for sample surveys}
\description{
Compute means, variances, quantiles, and cross-tabulations for data from complex surveys.
}
\usage{
svyquantile(x, design, quantiles, method = "linear", f = 1)
svymean(x, design, na.rm=FALSE) 
svrepmean(x, design, na.rm=FALSE, rho=NULL, return.replicates=FALSE) 
svyvar(x, design, na.rm=FALSE) 
svytotal(x, design, na.rm=FALSE) 
svreptotal(x, design, na.rm=FALSE, rho=NULL, return.replicates=FALSE) 
svytable(formula, design, Ntotal = design$fpc, round = FALSE)
svreptable(formula, design, Ntotal = sum(weights(design, "sampling"))), round = FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A formula, vector or matrix}
  \item{design}{\code{survey.design} object}
  \item{quantiles}{Quantiles to estimate}
  \item{method}{see \code{\link{approxfun}}}
  \item{f}{see \code{\link{approxfun}}}
  \item{na.rm}{Should missing values be removed?}
  \item{formula}{A one-sided formula specifying variables to be tabulated}
  \item{Ntotal}{A population total or set of population stratum totals
    to normalise to.}
  \item{round}{Should the table entries be rounded to the nearest
    integer?}
  \item{rho}{parameter for Fay's variance estimator in a BRR design}
  \item{return.replicates}{Return the replicate means?}
}
\details{
  
These functions perform weighted estimation, with each observation being
weighted by the inverse of its sampling probability.  The \code{svymean}
and \code{svyvar} functions also give precision estimates that
incorporate the effects of stratification and clustering. The first
four functions are similar to the standard functions whose names do not
begin with \code{svy}.

The \code{svytotal} and \code{svreptotal} functions estimate a
population total.  Use \code{predict} on \code{\link{svyratio}},
\code{\link{svrepratio}}, \code{\link{svyglm}}, \code{\link{svrepglm}}
to get ratio or regression estimates of totals.

The \code{svytable} and \code{svreptable} function computes a weighted
crosstabulation.  If the sampling probabilities supplied to
\code{svydesign} were actual probabilities (rather than relative
probabilities) this estimates a full population crosstabulation.
Otherwise it estimates only relative proportions and should be
normalised to some convenient total such as 100 or 1.0 by specifying the
\code{Ntotal} argument.

The \code{Ntotal} argument can be either a single number or a data frame
whose first column is the sampling strata and second column the
population size in each stratum.  In this second case the
\code{svytable} command performs `post-stratification': tabulating
and scaling to the population within strata and then adding up the
strata.

As with other \code{xtabs} objects, the output of \code{svytable} can be
processed by \code{ftable} for more attractive display.
}
\value{
  The first three functions return vectors, the last returns an \code{\link{xtabs}} object. 
}

\author{Thomas Lumley}

\seealso{ \code{\link{svydesign}}, \code{\link{svyCprod}}, \code{\link{mean}},\code{\link{var}}, \code{\link{quantile}}, \code{\link{xtabs}} }

\examples{

  data(api)
  ## population
  mean(apipop$api00)
  quantile(apipop$api00,c(.25,.5,.75))
  var(apipop$api00)
  sum(apipop$enroll)
  sum(apipop$api.stu)/sum(apipop$enroll)
  xtabs(~sch.wide+comp.imp, data=apipop)

  ## one-stage cluster sample
  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
  summary(dclus1)
  svymean(~api00, dclus1)
  svyquantile(~api00, dclus1, c(.25,.5,.75))
  svyvar(~api00, dclus1)
  svytotal(~enroll, dclus1)
  svyratio(~api.stu, ~enroll, dclus1)
  svytable(~sch.wide+comp.imp, dclus1)


  #stratified sample
  dstrat<-svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  summary(dstrat)
  svymean(~api00, dstrat)
  svyquantile(~api00, dstrat, c(.25,.5,.75))
  svyvar(~api00, dstrat)
  svytotal(~enroll, dstrat)
  svyratio(~api.stu, ~enroll, dstrat)
  svytable(~sch.wide+comp.imp, dstrat)

  data(scd)
  repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
              c(0,1,0,1,1,0))
  scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)
  svrepmean(~arrests+alive, design=scdrep)

 }
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
