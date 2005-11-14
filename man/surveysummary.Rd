\name{surveysummary}
\alias{svymean}
\alias{svymean.survey.design}
\alias{svymean.survey.design2}
\alias{svymean.svyrep.design}
\alias{svrepmean}
\alias{svytotal}
\alias{svytotal.survey.design}
\alias{svytotal.survey.design2}
\alias{svytotal.svyrep.design}
\alias{print.svystat}
\alias{print.svrepstat}
\alias{svreptotal}
\alias{svyvar}
\alias{svyvar.survey.design}
\alias{svyvar.svyrep.design}
\alias{coef.svystat}
\alias{vcov.svystat}
\alias{coef.svrepstat}
\alias{vcov.svrepstat}
\alias{as.data.frame.svystat}
\alias{as.data.frame.svrepstat}
\alias{cv.svyratio}
\alias{cv.svrepratio}
\alias{cv.svrepstat}
\alias{cv.svystat}
\alias{cv.default}
\alias{cv}
\alias{deff}
\alias{deff.default}
\alias{svrepvar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summary statistics for sample surveys}
\description{
Compute means, variances, ratios and totals for data from complex surveys.
}
\usage{
\method{svymean}{survey.design}(x, design, na.rm=FALSE,deff=FALSE,...) 
\method{svymean}{svyrep.design}(x, design, na.rm=FALSE, rho=NULL,
  return.replicates=FALSE, deff=FALSE,...) 
\method{svyvar}{survey.design}(x, design, na.rm=FALSE,...) 
\method{svyvar}{svyrep.design}(x, design, na.rm=FALSE, rho=NULL,
   return.replicates=FALSE,...) 
\method{svytotal}{survey.design}(x, design, na.rm=FALSE,deff=FALSE,...) 
\method{svytotal}{svyrep.design}(x, design, na.rm=FALSE, rho=NULL,
   return.replicates=FALSE, deff=FALSE,...)
\method{coef}{svystat}(object,...)
\method{coef}{svrepstat}(object,...)
\method{vcov}{svystat}(object,...)
\method{vcov}{svrepstat}(object,...)
cv(object,...)
deff(object, quietly=FALSE,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A formula, vector or matrix}
  \item{design}{\code{survey.design} or \code{svyrep.design} object}
  \item{na.rm}{Should cases with missing values be dropped?}
  \item{rho}{parameter for Fay's variance estimator in a BRR design}
  \item{return.replicates}{Return the replicate means?}
  \item{deff}{Return the design effect (see below)}
  \item{object}{The result of one of the other survey summary functions}
  \item{quietly}{Don't warn when there is no design effect computed}
  \item{...}{additional arguments to \code{cv} methods,not currently used}
}
\details{
  
These functions perform weighted estimation, with each observation being
weighted by the inverse of its sampling probability.  Except for the
table functions, these also give precision estimates that incorporate
the effects of stratification and clustering.

Factor variables are converted to sets of indicator variables for each
category in computing means and totals. Combining this with the
\code{\link{interaction}} function, allows crosstabulations. See
\code{\link{ftable.svystat}} for formatting the output.

With \code{na.rm=TRUE}, all cases with missing data are removed. With
\code{na.rm=FALSE} cases with missing data are not removed and so will
produce missing results.  When using replicate weights and
\code{na.rm=FALSE} it may be useful to set
\code{options(na.action="na.pass")}, otherwise all replicates with any
missing results will be discarded.

The \code{svytotal} and \code{svreptotal} functions estimate a
population total.  Use \code{predict} on \code{\link{svyratio}},
\code{\link{svrepratio}}, \code{\link{svyglm}}, \code{\link{svrepglm}}
to get ratio or regression estimates of totals.

The design effect compares the variance of a mean or total to the
variance from a study of the same size using simple random sampling
without replacement. Note that the design effect will be incorrect if
the weights have been rescaled so that they are not reciprocals of
sampling probabilities.  To obtain an estimate of the design effect
comparing to simple random sampling with replacement, which does not
have this requirement, use \code{deff="replace"}. This with-replacement
design effect is the square of Kish's "deft".


The \code{cv} function computes the coefficient of variation of a
statistic such as ratio, mean or total. The default method is for any
object with methods for \code{\link{SE}} and \code{coef}.

}
\value{
  Objects of class \code{"svystat"} or \code{"svrepstat"},
    which are vectors with a \code{"var"} attribute giving the variance
    and a \code{"statistic"} attribute giving the name of the statistic.
    
}

\author{Thomas Lumley}

\seealso{ \code{\link{svydesign}}, \code{\link{as.svrepdesign}},
  \code{\link{svrepdesign}},  \code{\link{svyquantile}} ,  \code{\link{ftable.svystat}}}

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
  svymean(~api00, dclus1, deff=TRUE)
  svymean(~factor(stype),dclus1)
  svymean(~interaction(stype, comp.imp), dclus1)
  svyquantile(~api00, dclus1, c(.25,.5,.75))
  svyvar(~api00, dclus1)
  svytotal(~enroll, dclus1, deff=TRUE)
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
  svymean(~api00, jkstrat)
  svymean(~factor(stype),jkstrat)
  svyvar(~api00,jkstrat)
  svyquantile(~api00, jkstrat, c(.25,.5,.75))
  svytotal(~enroll, jkstrat)
  svyratio(~api.stu, ~enroll, jkstrat)

  # coefficients of variation
  cv(svytotal(~enroll,dstrat))
  cv(svyratio(~api.stu, ~enroll, jkstrat))

  # extracting statistic and variance
  coef(svytotal(~enroll,dstrat))
  vcov(svymean(~api00+api99,jkstrat))

  # Design effect
  svymean(~api00, dstrat, deff=TRUE)
  svymean(~api00, dstrat, deff="replace")
  svymean(~api00, jkstrat, deff=TRUE)
  svymean(~api00, jkstrat, deff="replace")
 (a<-svytotal(~enroll, dclus1, deff=TRUE))
  deff(a)

  # BRR method
  data(scd)
  repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
              c(0,1,0,1,1,0))
  scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)
  svymean(~arrests+alive, design=scdrep)

 }
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
