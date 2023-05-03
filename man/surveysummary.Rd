\name{surveysummary}
\alias{svymean}
\alias{svymean.survey.design}
\alias{svymean.survey.design2}
\alias{svymean.svyrep.design}
\alias{svymean.twophase}
\alias{svytotal}
\alias{svytotal.twophase}
\alias{svytotal.survey.design}
\alias{svytotal.survey.design2}
\alias{svytotal.svyrep.design}
\alias{svyvar}
\alias{svyvar.survey.design}
\alias{svyvar.svyrep.design}
\alias{coef.svystat}
\alias{vcov.svystat}
\alias{coef.svrepstat}
\alias{vcov.svrepstat}
\alias{cv.svyratio}
\alias{cv.svrepratio}
\alias{cv.svrepstat}
\alias{cv.svystat}
\alias{cv.default}
\alias{cv}
\alias{deff}
\alias{deff.default}
\alias{confint.svystat}
\alias{confint.svrepstat}
\alias{make.formula}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summary statistics for sample surveys}
\description{
Compute means, variances, ratios and totals for data from complex surveys.
}
\usage{
\method{svymean}{survey.design}(x, design, na.rm=FALSE,deff=FALSE,influence=FALSE,...) 
\method{svymean}{survey.design2}(x, design, na.rm=FALSE,deff=FALSE,influence=FALSE,...) 
\method{svymean}{twophase}(x, design, na.rm=FALSE,deff=FALSE,...) 
\method{svymean}{svyrep.design}(x, design, na.rm=FALSE, rho=NULL,
  return.replicates=FALSE, deff=FALSE,...) 
\method{svyvar}{survey.design}(x, design, na.rm=FALSE,...) 
\method{svyvar}{svyrep.design}(x, design, na.rm=FALSE, rho=NULL,
   return.replicates=FALSE,...,estimate.only=FALSE) 
\method{svytotal}{survey.design}(x, design, na.rm=FALSE,deff=FALSE,influence=FALSE,...) 
\method{svytotal}{survey.design2}(x, design, na.rm=FALSE,deff=FALSE,influence=FALSE,...) 
\method{svytotal}{twophase}(x, design, na.rm=FALSE,deff=FALSE,...) 
\method{svytotal}{svyrep.design}(x, design, na.rm=FALSE, rho=NULL,
   return.replicates=FALSE, deff=FALSE,...)
\method{coef}{svystat}(object,...)
\method{coef}{svrepstat}(object,...)
\method{vcov}{svystat}(object,...)
\method{vcov}{svrepstat}(object,...)
\method{confint}{svystat}(object,  parm, level = 0.95,df =Inf,...)
\method{confint}{svrepstat}(object,  parm, level = 0.95,df =Inf,...)
cv(object,...)
deff(object, quietly=FALSE,...)
make.formula(names)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A formula, vector or matrix}
  \item{design}{\code{survey.design} or \code{svyrep.design} object}
  \item{na.rm}{Should cases with missing values be dropped?}
  \item{influence}{Should a matrix of influence functions be returned
    (primarily to support \code{\link{svyby}})}
  \item{rho}{parameter for Fay's variance estimator in a BRR design}
  \item{return.replicates}{Return the replicate means/totals?}
  \item{deff}{Return the design effect (see below)}
  \item{object}{The result of one of the other survey summary functions}
  \item{quietly}{Don't warn when there is no design effect computed}
  \item{estimate.only}{Don't compute standard errors (useful when
    \code{svyvar} is used to estimate the design effect)}
  \item{parm}{a specification of which parameters are to be given
    confidence intervals, either a vector of numbers or a vector of
    names. If missing, all parameters are considered.}
  \item{level}{the confidence level required.}
  \item{df}{degrees of freedom for t-distribution in confidence
    interval, use \code{degf(design)} for number of PSUs minus number of
  strata}
  \item{...}{additional arguments to  methods,not currently
    used}
  \item{names}{vector of character strings}
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
population total.  Use \code{predict} on \code{\link{svyratio}} and
\code{\link{svyglm}}, to get ratio or regression estimates of totals.

\code{svyvar} estimates the population variance. The object returned
includes the full matrix of estimated population variances and
covariances, but by default only the diagonal elements are printed. To
display the whole matrix use \code{as.matrix(v)} or \code{print(v,
covariance=TRUE)}.

The design effect compares the variance of a mean or total to the
variance from a study of the same size using simple random sampling
without replacement. Note that the design effect will be incorrect if
the weights have been rescaled so that they are not reciprocals of
sampling probabilities.  To obtain an estimate of the design effect
comparing to simple random sampling with replacement, which does not
have this requirement, use \code{deff="replace"}. This with-replacement
design effect is the square of Kish's "deft".

The design effect for a subset of a design conditions on the size of
the subset. That is, it compares the variance of the estimate to the
variance of an estimate based on a simple random sample of the same
size as the subset, taken from the subpopulation. So, for example,
under stratified random sampling the design effect in a subset
consisting of a single stratum will be 1.0.


The \code{cv} function computes the coefficient of variation of a
statistic such as ratio, mean or total. The default method is for any
object with methods for \code{\link{SE}} and \code{coef}.

\code{make.formula} makes a formula from a vector of names.  This is
useful because formulas as the best way to specify variables to the
survey functions.

}
\value{
  Objects of class \code{"svystat"} or \code{"svrepstat"},
    which are vectors with a \code{"var"} attribute giving the variance
    and a \code{"statistic"} attribute giving the name of the
    statistic.

    These objects have methods for \code{vcov}, \code{SE}, \code{coef},
    \code{confint}, \code{svycontrast}.
    
}

\author{Thomas Lumley}

\seealso{ \code{\link{svydesign}}, \code{\link{as.svrepdesign}},
  \code{\link{svrepdesign}} for constructing design objects.

  \code{\link{degf}} to extract degrees of freedom from a design.

  \code{\link{svyquantile}}  for quantiles
  
  \code{\link{ftable.svystat}} for more attractive tables

  \code{\link{svyciprop}} for more accurate confidence intervals for
  proportions near 0 or 1.

  \code{\link{svyttest}} for comparing two means.

  \code{\link{svycontrast}} for linear and nonlinear functions of estimates.
}

\examples{

  data(api)

  ## one-stage cluster sample
  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

  svymean(~api00, dclus1, deff=TRUE)
  svymean(~factor(stype),dclus1)
  svymean(~interaction(stype, comp.imp), dclus1)
  svyquantile(~api00, dclus1, c(.25,.5,.75))
  svytotal(~enroll, dclus1, deff=TRUE)
  svyratio(~api.stu, ~enroll, dclus1)

  v<-svyvar(~api00+api99, dclus1)
  v
  print(v, cov=TRUE)
  as.matrix(v)


  # replicate weights - jackknife (this is slower)
  dstrat<-svydesign(id=~1,strata=~stype, weights=~pw,
        data=apistrat, fpc=~fpc)
  jkstrat<-as.svrepdesign(dstrat)

  svymean(~api00, jkstrat)
  svymean(~factor(stype),jkstrat)
  svyvar(~api00+api99,jkstrat)

  svyquantile(~api00, jkstrat, c(.25,.5,.75))
  svytotal(~enroll, jkstrat)
  svyratio(~api.stu, ~enroll, jkstrat)

  # coefficients of variation
  cv(svytotal(~enroll,dstrat))
  cv(svyratio(~api.stu, ~enroll, jkstrat))

  # extracting information from the results
  coef(svytotal(~enroll,dstrat))
  vcov(svymean(~api00+api99,jkstrat))
  SE(svymean(~enroll, dstrat))
  confint(svymean(~api00+api00, dclus1))
  confint(svymean(~api00+api00, dclus1), df=degf(dclus1))

  # Design effect
  svymean(~api00, dstrat, deff=TRUE)
  svymean(~api00, dstrat, deff="replace")
  svymean(~api00, jkstrat, deff=TRUE)
  svymean(~api00, jkstrat, deff="replace")
 (a<-svytotal(~enroll, dclus1, deff=TRUE))
  deff(a)

## weights that are *already* calibrated to population size
sum(weights(dclus1))
nrow(apipop)
cdclus1<- svydesign(id=~dnum, weights=~pw, data=apiclus1,
fpc=~fpc,calibrate.formula=~1)
SE(svymean(~enroll, dclus1))
## not equal to SE(mean)
SE(svytotal(~enroll, dclus1))/nrow(apipop)
## equal to SE(mean)
SE(svytotal(~enroll, cdclus1))/nrow(apipop)

 }

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
