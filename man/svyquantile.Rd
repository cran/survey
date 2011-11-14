\name{svyquantile}
\alias{svyquantile}
\alias{print.svyquantile}
\alias{SE.svyquantile}
\alias{svyquantile.survey.design}
\alias{svyquantile.svyrep.design}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Quantiles for sample surveys}
\description{
Compute quantiles for data from complex surveys.
}
\usage{
\method{svyquantile}{survey.design}(x, design, quantiles, alpha=0.05,
   ci=FALSE, method = "linear", f = 1,
   interval.type=c("Wald","score","betaWald"), na.rm=FALSE,se=ci,
   ties=c("discrete","rounded"), df=Inf,...)
\method{svyquantile}{svyrep.design}(x, design, quantiles,
   method ="linear", interval.type=c("probability","quantile"), f = 1,
   return.replicates=FALSE, ties=c("discrete","rounded"),na.rm=FALSE,...)
\method{SE}{svyquantile}(object,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A formula, vector or matrix}
  \item{design}{\code{survey.design} or \code{svyrep.design} object}
  \item{quantiles}{Quantiles to estimate}
  \item{method}{see \code{\link{approxfun}}}
  \item{f}{see \code{\link{approxfun}}}
  \item{ci}{Compute a confidence interval? (relatively slow; needed for \code{\link{svyby}})}
  \item{se}{Compute standard errors from the confidence interval length?}
  \item{alpha}{Level for confidence interval}
  \item{interval.type}{See Details below}
  \item{ties}{See Details below}
  \item{df}{Degrees of freedom for a t-distribution. \code{Inf} requests a Normal distribution,
   \code{NULL} uses \code{\link{degf}}. Not relevant for \code{type="betaWald"}}
  \item{return.replicates}{Return the replicate means?}
  \item{na.rm}{Remove \code{NA}s?}
  \item{...}{arguments for future expansion}
  \item{object}{Object returned by \code{svyquantile.survey.design}}
}
\details{
   The definition of the CDF and thus of the quantiles is ambiguous in
  the presence of ties.  With \code{ties="discrete"} the data are
  treated as genuinely discrete, so the CDF has vertical steps at tied
  observations. With \code{ties="rounded"} all the weights for tied
  observations are summed and the CDF interpolates linearly between
  distinct observed values, and so is a continuous function.  Combining
  \code{interval.type="betaWald"} and \code{ties="discrete"} is (close
  to) the proposal of Shah and Vaish(2006) used in some versions of SUDAAN.
  
  Interval estimation for quantiles is complicated, because the
  influence function is not continuous.  Linearisation cannot be used
  directly, and computing the variance of replicates is valid only for
  some designs (eg BRR, but not jackknife). The \code{interval.type}
  option controls how the intervals are computed.
  
  For \code{survey.design} objects the default is
  \code{interval.type="Wald"}. A 95\% Wald confidence interval is
  constructed for the proportion below the estimated quantile. The
  inverse of the estimated CDF is used to map this to a confidence
  interval for the quantile. This is the method of Woodruff
  (1952). For \code{"betaWald"} the same procedure is used, but the
  confidence interval for the proportion is computed using the exact
  binomial cdf with an effective sample size proposed by Korn &
  Graubard (1998).


  If \code{interval.type="score"} we use a method described by Binder
  (1991) and due originally to Francisco and Fuller (1986), which
  corresponds to inverting a robust score test.  At the upper and lower
  limits of the confidence interval, a test of the null hypothesis that
  the cumulative distribution function is equal to the target quantile
  just rejects.  This was the default before version 2.9. It is much
  slower than \code{"Wald"}, and Dorfman & Valliant (1993) suggest it is
  not any more accurate.
  
  Standard errors are computed from these confidence intervals by
  dividing the confidence interval length by \code{2*qnorm(alpha/2)}.
  
  For replicate-weight designs, ordinary replication-based standard errors
  are valid for BRR and Fay's method, and for some bootstrap-based
  designs, but not for jackknife-based designs.
  \code{interval.type="quantile"} gives these replication-based
  standard errors.  The default, \code{interval.type="probability"}
  computes confidence on the probability scale and then transforms
  back to quantiles, the equivalent of \code{interval.type="Wald"} for
  \code{survey.design} objects (with \code{alpha=0.05}).

  There is a \code{confint} method for \code{svyquantile} objects; it
  simply extracts the pre-computed confidence interval.
  
}
\value{
  returns a list whose first component is the quantiles and second
  component is the confidence intervals. For replicate weight designs,
  returns an object of class \code{svyrepstat}.
}


\author{Thomas Lumley}


\seealso{
  \code{\link{svykm}} for quantiles of survival curves

  \code{\link{svyciprop}} for confidence intervals on proportions.
}
\references{
  Binder DA (1991) Use of estimating functions for interval estimation
  from complex surveys. \emph{Proceedings of the ASA Survey Research
    Methods Section}  1991: 34-42

  Dorfman A, Valliant R (1993) Quantile variance estimators in complex
  surveys. Proceedings of the ASA Survey Research Methods Section. 1993: 866-871

  Korn EL, Graubard BI. (1998) Confidence Intervals For Proportions With
  Small Expected Number of Positive Counts Estimated From Survey
  Data. Survey Methodology 23:193-201.
  
  Francisco CA, Fuller WA (1986) Estimation of the distribution
  function with a complex survey. Technical Report, Iowa State
  University.
  
  Shao J, Tu D (1995) \emph{The Jackknife and Bootstrap}. Springer.

  Shah BV, Vaish AK (2006) Confidence Intervals for Quantile Estimation
  from Complex Survey Data. Proceedings of the Section on Survey
  Research Methods. 
  
  Woodruff RS (1952) Confidence intervals for medians and other
  position measures. JASA 57, 622-627.  
}
\examples{

  data(api)
  ## population
  quantile(apipop$api00,c(.25,.5,.75))

  ## one-stage cluster sample
  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
  svyquantile(~api00, dclus1, c(.25,.5,.75),ci=TRUE)
  svyquantile(~api00, dclus1, c(.25,.5,.75),ci=TRUE,interval.type="betaWald")
  svyquantile(~api00, dclus1, c(.25,.5,.75),ci=TRUE,df=NULL)

  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
  (qapi<-svyquantile(~api00, dclus1, c(.25,.5,.75),ci=TRUE, interval.type="score"))
  SE(qapi)

  #stratified sample
  dstrat<-svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  svyquantile(~api00, dstrat, c(.25,.5,.75),ci=TRUE)

  #stratified sample, replicate weights
  # interval="probability" is necessary for jackknife weights
  rstrat<-as.svrepdesign(dstrat)
  svyquantile(~api00, rstrat, c(.25,.5,.75), interval.type="probability")


  # BRR method
  data(scd)
  repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
              c(0,1,0,1,1,0))
  scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)
  svyquantile(~arrests+alive, design=scdrep, quantile=0.5, interval.type="quantile")

 }
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
