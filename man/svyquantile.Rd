\name{newsvyquantile}
\alias{svyquantile}
\alias{svyquantile.survey.design}
\alias{svyquantile.svyrep.design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantiles under complex sampling.
}
\description{
Estimates quantiles and confidence intervals for them. This function was
completely re-written for version 4.1 of the survey package, and has a
wider range of ways to define the quantile. See the vignette for a list of them. 
}
\usage{
svyquantile(x, design, quantiles, ...)
\method{svyquantile}{survey.design}(x, design, quantiles, alpha = 0.05,
                      interval.type = c("mean", "beta","xlogit", "asin","score"),
                      na.rm = FALSE,  ci=TRUE, se = ci,
                      qrule=c("math","school","shahvaish","hf1","hf2","hf3",
		      "hf4","hf5","hf6","hf7","hf8","hf9"),
                      df = NULL, ...)
\method{svyquantile}{svyrep.design}(x, design, quantiles, alpha = 0.05,
                      interval.type = c("mean", "beta","xlogit", "asin","quantile"),
                      na.rm = FALSE, ci = TRUE, se=ci,
                      qrule=c("math","school","shahvaish","hf1","hf2","hf3",
		      "hf4","hf5","hf6","hf7","hf8","hf9"),
                      df = NULL, return.replicates=FALSE,...)		      
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A one-sided formula describing variables to be used
}
  \item{design}{
Design object
}
  \item{quantiles}{
Numeric vector specifying which quantiles are requested
}
\item{alpha}{Specified confidence interval coverage}
\item{interval.type}{See Details below}
\item{na.rm}{Remove missing values?}
\item{ci,se}{Return an estimated confidence interval and standard error?}
\item{qrule}{Rule for defining the quantiles: either a character string specifying one of the built-in rules, or a function}
\item{df}{Degrees of freedom for confidence interval estimation: \code{NULL} specifies \code{degf(design)}}
\item{return.replicates}{Return replicate estimates of the quantile (only for \code{interval.type="quantile"})}
  \item{\dots}{
For future expansion
}
}
\details{
The \code{p}th quantile is defined as the value where the estimated cumulative
distribution function is equal to \code{p}. As with quantiles in
unweighted data, this definition only pins down the quantile to an
interval between two observations, and a rule is needed to interpolate.
The default is the mathematical definition, the lower end of the
quantile interval; \code{qrule="school"} uses the midpoint of the
quantile interval; \code{"hf1"} to \code{"hf9"} are weighted analogues of
\code{type=1} to \code{9} in \code{\link{quantile}}. See the vignette
"Quantile rules" for details and for how to write your own. 
  
By default, confidence intervals are estimated using Woodruff's (1952) method,
which involves computing the quantile, estimating a confidence interval
for the proportion of observations below the quantile, and then
transforming that interval using the estimated CDF.  In that context,
the \code{interval.type} argument specifies how the confidence interval
for the proportion is computed, matching \code{\link{svyciprop}}. In
contrast to \code{oldsvyquantile}, \code{NaN} is returned if a confidence
interval endpoint on the probability scale falls outside \code{[0,1]}.

There are two exceptions. For \code{svydesign} objects,
\code{interval.type="score"} asks for the Francisco & Fuller confidence
interval based on inverting a score test. According to Dorfmann &
Valliant, this interval has inferior performance to the \code{"beta"}
and \code{"logit"} intervals; it is provided for compatibility.

For replicate-weight designs, \code{interval.type="quantile"} ask for an
interval based directly on the replicates of the quantile. This interval
is not valid for jackknife-type replicates, though it should perform well for
bootstrap-type replicates, BRR, and SDR.

The \code{df} argument specifies degrees of freedom for a t-distribution
approximation to distributions of means. The default is the design degrees of
freedom. Specify \code{df=Inf} to use a Normal distribution (eg, for compatibility).

When the standard error is requested, it is estimated by dividing the
confidence interval length by the number of standard errors in a t
confidence interval with the specified \code{alpha}. For example, with
\code{alpha=0.05} and \code{df=Inf} the standard error is estimated as the confidence
interval length divided by \code{2*1.96}.

}
\value{
An object of class \code{"newsvyquantile"}, except that with a
replicate-weights design and \code{interval.type="quantile"} and
\code{return.replicates=TRUE} it's an object of class \code{"svrepstat"}
}
\references{
  Dorfman A, Valliant R (1993) Quantile variance estimators in complex
  surveys. Proceedings of the ASA Survey Research Methods Section. 1993: 866-871
  
  Francisco CA, Fuller WA (1986) Estimation of the distribution
  function with a complex survey. Technical Report, Iowa State
  University.

  Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical packages,
  The American Statistician 50, 361-365.

  Shah BV, Vaish AK (2006) Confidence Intervals for Quantile Estimation
  from Complex Survey Data. Proceedings of the Section on Survey
  Research Methods. 
  
  Woodruff RS (1952) Confidence intervals for medians and other
  position measures. JASA 57, 622-627.
}



\seealso{
  \code{vignette("qrule", package = "survey")}
  \code{\link{oldsvyquantile}}	 
\code{\link{quantile}}
}
\examples{
data(api)
## population
quantile(apipop$api00,c(.25,.5,.75))

## one-stage cluster sample
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1<-as.svrepdesign(dclus1)
bclus1<-as.svrepdesign(dclus1,type="boot")


svyquantile(~api00, dclus1, c(.25,.5,.75))
svyquantile(~api00, dclus1, c(.25,.5,.75),interval.type="beta")

svyquantile(~api00, rclus1, c(.25,.5,.75))
svyquantile(~api00, rclus1, c(.25,.5,.75),interval.type="quantile")
svyquantile(~api00, bclus1, c(.25,.5,.75),interval.type="quantile")

svyquantile(~api00+ell, dclus1, c(.25,.5,.75), qrule="math")
svyquantile(~api00+ell, dclus1, c(.25,.5,.75), qrule="school")
svyquantile(~api00+ell, dclus1, c(.25,.5,.75), qrule="hf8")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}% use one of  RShowDoc("KEYWORDS")
\keyword{univar}% at least one, from doc/KEYWORDS
