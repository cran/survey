\name{svyquantile}
\alias{svyquantile}
\alias{print.svyquantile}
\alias{SE.svyquantile}
\alias{svyquantile.survey.design}
\alias{svyquantile.svyrep.design}
\alias{svrepquantile}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Quantiles for sample surveys}
\description{
Compute quantiles for data from complex surveys.
}
\usage{
\method{svyquantile}{survey.design}(x, design, quantiles, alpha=0.05, ci=FALSE, method = "linear", f = 1, interval.type=c("Wald","score"),...)
\method{svyquantile}{svyrep.design}(x, design, quantiles, method =
"linear", f = 1, return.replicates=FALSE,...)
\method{SE}{svyquantile}(object,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A formula, vector or matrix}
  \item{design}{\code{survey.design} or \code{svyrep.design} object}
  \item{quantiles}{Quantiles to estimate}
  \item{method}{see \code{\link{approxfun}}}
  \item{f}{see \code{\link{approxfun}}}
  \item{ci}{Compute a confidence interval (relatively slow)?}
  \item{alpha}{Level for confidence interval}
  \item{interval.type}{See Details below}
  \item{return.replicates}{Return the replicate means?}
  \item{...}{arguments for future expansion}
  \item{object}{Object returned by \code{svyquantile.survey.design}}
}
\details{
  Interval estimation for quantiles is complicated, because the
  influence function is not continuous.  Linearisation cannot be used,
  and only some replication weight designs give valid results.
  
  For \code{survey.design} objects the computation is controlled by
  \code{interval.type}. If \code{"Wald"}, a 95\% Wald confidence
  interval is constructed for the proportion below the estimated
  quantile. The inverse of the estimated CDF is used to map this to a
  confidence interval for the quantile.
  
  If \code{interval.type="score"} we use the method of Francisco-Fuller, which
  corresponds to inverting a robust score test.  At the upper and lower
  limits of the confidence interval, a test of the null hypothesis that
  the cumulative distribution function is equal to the target quantile
  just rejects.  This was the default before version 2.9 and should be
  more accurate, but is much slower.
  
  Standard errors are computed from these confidence intervals by
  dividing the confidence interval length by \code{2*qnorm(alpha/2)}.
  
  For replicate-weight designs, ordinary replication-based standard errors
  are computed.  These are typically not valid for the JK1 and JKn jackknife
  designs. They are valid for BRR and Fay's method, and for some
  bootstrap-based designs.

  
}
\value{
  returns a list whose first component is the quantiles and second
  component is the confidence intervals. For replicate weight designs,
  returns an object of class \code{svyrepstat}.
}

\author{Thomas Lumley}

\seealso{ \code{\link{svydesign}}, \code{\link{svymean}},
  \code{\link{as.svrepdesign}}, \code{\link{svrepdesign}}}
\references{
  Binder DA (1991) Use of estimating functions for interval estimation
  from complex surveys. \emph{Proceedings of the ASA Survey Research
    Methods Section}  1991: 34-42
  
  Shao J, Tu D (1995) \emph{The Jackknife and Bootstrap}. Springer.  
  }
\examples{

  data(api)
  ## population
  quantile(apipop$api00,c(.25,.5,.75))

  ## one-stage cluster sample
  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
  svyquantile(~api00, dclus1, c(.25,.5,.75),ci=TRUE)

  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
  (qapi<-svyquantile(~api00, dclus1, c(.25,.5,.75),ci=TRUE, interval.type="score"))
  SE(qapi)

  #stratified sample
  dstrat<-svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  svyquantile(~api00, dstrat, c(.25,.5,.75),ci=TRUE)


  # BRR method
  data(scd)
  repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
              c(0,1,0,1,1,0))
  scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)
  svyquantile(~arrests+alive, design=scdrep, quantile=0.5)

 }
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
