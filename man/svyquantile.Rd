\name{svyquantile}
\alias{svyquantile}
\alias{svrepquantile}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Quantiles for sample surveys}
\description{
Compute quantiles for data from complex surveys.
}
\usage{
svyquantile(x, design, quantiles, alpha=0.05, ci=FALSE,method = "linear", f = 1)
svrepquantile(x, design, quantiles, method = "linear", f = 1, return.replicates=FALSE)
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
  \item{return.replicates}{Return the replicate means?}
}
\details{
  Interval estimation for quantiles is complicated, because the
  influence function is not continuous.  Linearisation cannot be used,
  and only some replication weight designs give valid results.
  
  For \code{svyrepquantile} we use the method of Francisco-Fuller, which
  corresponds to inverting a robust score test.  At the upper and lower
  limits of the confidence interval, a test of the null hypothesis that
  the cumulative distribution function is equal to the target quantile
  just rejects.
  
  For \code{svrepquantile} ordinary replication-based standard errors
  are computed.  These are not valid for the JK1 and JKn jackknife
  designs. They are valid for BRR and Fay's method, and for some
  bootstrap-based designs.

  
}
\value{
  \code{svyquantile} returns a list whose first component is the
  quantiles and second component is the confidence
  intervals. \code{svrepquantile} returns an object of class \code{svyrepstat}.
}

\author{Thomas Lumley}

\seealso{ \code{\link{svydesign}}, \code{\link{svymean}},
  \code{\link{as.svrepdesign}}, \code{\link{svrepdesign}}}
\references{
  Binder DA (1991) Use of estimating functions for interval estimation
  from complex surveys. \emph{Journal of Official Statistics}  1991: 34-42
  Shao J, Tu D (1995) \emph{The Jackknife and Bootstrap}. Springer.  
  }
\examples{

  data(api)
  ## population
  quantile(apipop$api00,c(.25,.5,.75))

  ## one-stage cluster sample
  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
  svyquantile(~api00, dclus1, c(.25,.5,.75),ci=TRUE)

  #stratified sample
  dstrat<-svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  svyquantile(~api00, dstrat, c(.25,.5,.75),ci=TRUE)


  # BRR method
  data(scd)
  repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
              c(0,1,0,1,1,0))
  scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)
  svrepquantile(~arrests+alive, design=scdrep, quantile=0.5)

 }
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
