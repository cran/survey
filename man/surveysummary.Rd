\name{surveysummary}
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
Compute means, variances, ratios and totals for data from complex surveys.
}
\usage{
svymean(x, design, na.rm=FALSE) 
svrepmean(x, design, na.rm=FALSE, rho=NULL, return.replicates=FALSE) 
svyvar(x, design, na.rm=FALSE) 
svytotal(x, design, na.rm=FALSE) 
svreptotal(x, design, na.rm=FALSE, rho=NULL, return.replicates=FALSE) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A formula, vector or matrix}
  \item{design}{\code{survey.design} or \code{svyrep.design} object}
  \item{na.rm}{Should missing values be removed?}
  \item{rho}{parameter for Fay's variance estimator in a BRR design}
  \item{return.replicates}{Return the replicate means?}
}
\details{
  
These functions perform weighted estimation, with each observation being
weighted by the inverse of its sampling probability.  Except for the
table functions, these also give precision estimates that incorporate
the effects of stratification and clustering.

The \code{svytotal} and \code{svreptotal} functions estimate a
population total.  Use \code{predict} on \code{\link{svyratio}},
\code{\link{svrepratio}}, \code{\link{svyglm}}, \code{\link{svrepglm}}
to get ratio or regression estimates of totals.

}
\value{
  Objects of class \code{"svystat"} or \code{"svrepstat"},
    which are vectors with a \code{"var"} attribute giving the variance
    and a \code{"statistic"} attribute giving the name of the statistic.
    
}

\author{Thomas Lumley}

\seealso{ \code{\link{svydesign}}, \code{\link{as.svrepdesign}},
  \code{\link{svrepdesign}}, \code{\link{svyCprod}}, \code{\link{mean}},
  \code{\link{var}}, \code{\link{svyquantile}} }

\examples{

  data(api)
  ## population
  mean(apipop$api00)
  quantile(apipop$api00,c(.25,.5,.75))
  var(apipop$api00)
  sum(apipop$enroll)
  sum(apipop$api.stu)/sum(apipop$enroll)

  ## one-stage cluster sample
  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
  summary(dclus1)
  svymean(~api00, dclus1)
  svyquantile(~api00, dclus1, c(.25,.5,.75))
  svyvar(~api00, dclus1)
  svytotal(~enroll, dclus1)
  svyratio(~api.stu, ~enroll, dclus1)

  #stratified sample
  dstrat<-svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  summary(dstrat)
  svymean(~api00, dstrat)
  svyquantile(~api00, dstrat, c(.25,.5,.75))
  svyvar(~api00, dstrat)
  svytotal(~enroll, dstrat)
  svyratio(~api.stu, ~enroll, dstrat)

  # replicate weights - jackknife (this is slow)
  jkstrat<-as.svrepdesign(dstrat)
  summary(jkstrat)
  svrepmean(~api00, jkstrat)
  svrepquantile(~api00, jkstrat, c(.25,.5,.75))
  svreptotal(~enroll, jkstrat)
  svrepratio(~api.stu, ~enroll, jkstrat)

  # BRR method
  data(scd)
  repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
              c(0,1,0,1,1,0))
  scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)
  svrepmean(~arrests+alive, design=scdrep)

 }
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
