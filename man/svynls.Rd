\name{svynls}
\alias{svynls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability-weighted nonlinear least squares
}
\description{
Fits a nonlinear model by probability-weighted least squares.  Uses
\code{nls} to do the fitting, but estimates design-based standard errors with either
linearisation or replicate weights. See \code{\link{nls}} for 
documentation of model specification and fitting.
}
\usage{
svynls(formula, design, start, weights=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    Nonlinear model specified as a formula; see \code{\link{nls}}
}
  \item{design}{
    Survey design object
  }
\item{start}{starting values, passed to \code{\link{nls}}}

  \item{weights}{
    Non-sampling weights, eg precision weights to give more efficient estimation in the presence of heteroscedasticity.
 }
  \item{\dots}{
    Other arguments to \code{nls} (especially, \code{start}). Also
    supports \code{return.replicates} for replicate-weight designs and
    \code{influence} for other designs. 
  }
}

\value{
  Object of class \code{svynls}. The fitted \code{nls} object is
  included as the \code{fit} element.
}

\seealso{
\code{\link{svymle}} for maximum likelihood with linear predictors on
one or more parameters
}
\examples{
set.seed(2020-4-3)
x<-rep(seq(0,50,1),10)
y<-((runif(1,10,20)*x)/(runif(1,0,10)+x))+rnorm(510,0,1)

pop_model<-nls(y~a*x/(b+x), start=c(a=15,b=5))

df<-data.frame(x=x,y=y)
df$p<-ifelse((y-fitted(pop_model))*(x-mean(x))>0, .4,.1)

df$strata<-ifelse(df$p==.4,"a","b")

in_sample<-stratsample(df$strata, round(table(df$strat)*c(0.4,0.1)))

sdf<-df[in_sample,]
des<-svydesign(id=~1, strata=~strata, prob=~p, data=sdf)
pop_model
(biased_sample<-nls(y~a*x/(b+x),data=sdf, start=c(a=15,b=5)))
(corrected <- svynls(y~a*x/(b+x), design=des, start=c(a=15,b=5)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey }% use one of  RShowDoc("KEYWORDS")
