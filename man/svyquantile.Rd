\name{svyquantile}
\alias{svyquantile}
\alias{svytable}
\alias{svymean}
\alias{svyvar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summary statistics for sample surveys}
\description{
Compute means, variances, quantiles, and cross-tabulations for data from complex surveys.
}
\usage{
svyquantile(x, design, quantiles, method = "linear", f = 1)
svymean(x, design, na.rm=FALSE) 
svyvar(x, design, na.rm=FALSE) 
svytable(formula, design, Ntotal = design$fpc, round = FALSE) 
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
  \item{round}{Should the table entries be rounded to the nearest integer?}
}
\details{
  
These functions perform weighted estimation, with each observation being
weighted by the inverse of its sampling probability.  The \code{svymean}
and \code{svyvar} functions also give precision estimates that
incorporate the effects of stratification and clustering. The first
three functions are similar to the standard functions whose names do not
begin with \code{svy}.

The \code{svytable} function computes a weighted crosstabulation.  If
the sampling probabilities supplied to \code{svydesign} were actual
probabilities (rather than relative probabilities) this estimates a full
population crosstabulation.  Otherwise it estimates only relative
proportions and should be normalised to some convenient total such as
100 or 1.0 by specifying the \code{Ntotal} argument.

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
\references{ ~put references to the literature/web site here ~ }
\author{Thomas Lumley}

\seealso{ \code{\link{svydesign}}, \code{\link{svyCprod}}, \code{\link{mean}},\code{\link{var}}, \code{\link{quantile}}, \code{\link{xtabs}} }

\examples{
 #population
  df<-data.frame(x=rnorm(1000),z=rep(0:4,200))
  df$y<-with(df, 3+3*x*z)
  #sampling fraction
  df$p<-with(df, exp(x)/(1+exp(x)))
  #sample
  xi<-rbinom(1000,1,df$p)
  sdf<-df[xi==1,]
  
  #survey design object: independent sampling, 
  dxi<-svydesign(~0,~p,data=sdf) 
  dxi

  mean(df$x)		#right
  mean(sdf$x)		#wrong
  svymean(sdf$x,dxi)	#right

  var(df$x)		#right
  var(sdf$x)		#wrong
  svyvar(sdf$x,dxi)		#right

  quantile(df$x,c(0.025,0.5,0.975))  #right
  quantile(sdf$x,c(0.025,0.5,0.975))  #wrong
  svyquantile(sdf$x,design=dxi,c(0.025,0.5,0.975))  #right

  table(sdf$z)  # sample table
  svytable(~z, dxi, round=TRUE) # estimated population table

 }
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{survey}% __ONLY ONE__ keyword per line
