\name{svyglm}
\alias{svyglm}
\alias{summary.svyglm}
\alias{print.summary.svyglm}
\alias{vcov.svyglm}
\alias{residuals.svyglm}
\alias{coef.svyglm}
\alias{extractAIC.svyglm}
\alias{logLik.svyglm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Survey-weighted generalised linear models.}
\description{
Fit a generalised linear model to data from a complex survey design, with inverse-probability weighting and with standard errors corrected for cluster sampling.
}
\usage{
svyglm(formula, design, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula}
  \item{design}{ \code{survey.design} object. Must contain all variables in the formula}
  \item{\dots}{Other arguments passed to \code{glm} }
}
\details{
Standard errors for cluster-sampled designs are computed using a linearisation estimate (equivalent to the Huber/White sandwich formula in GEEs). 
This does not incorporate the finite-population correction, making standard errors potentially conservative.
}
\value{
  An object of class \code{svyglm}.
}
\references{ ~put references to the literature/web site here ~ }
\author{Thomas Lumley}


\seealso{ \code{\link{svydesign}}, \code{\link{glm}}}

\examples{
  ## Independent sampling
  df<-data.frame(x=rnorm(1000),z=rep(0:4,200))
  df$y<-with(df, 3+3*x*z)
  df$p<-with(df, exp(x)/(1+exp(x)))
  mpop<-lm(y~x+z, data=df) 
  xi<-rbinom(1000,1,df$p)
  sdf<-df[xi==1,]

  dxi<-svydesign(~0,~p,data=sdf)
  m<-svyglm(y~x+z,data=sdf,family="gaussian",design=dxi)
  m1<-lm(y~x+z, data=sdf)

  summary(m1) ##wrong
  summary(m)  ##right
  summary(mpop) ##whole population

  ##cluster sampling
  df$id<-rep(1:250,each=4)
  df$clustp<-by(df,list(df$id),function(d) min(exp(d$x*d$z)/(1+exp(d$x*d$z))))[df$id]
  mpop<-lm(y~x+z, data=df) 

  xi<-rbinom(250,1,df$clustp[4*(1:250)])
  sdf<-df[xi[df$id]==1,]

  dxi<-svydesign(~id,~clustp,data=sdf)
  m<-svyglm(y~x+z,data=sdf,family="gaussian",design=dxi)
  m1<-lm(y~x+z,data=sdf)
 
  summary(m1) ##wrong
  summary(m)  ##right
  summary(mpop) ##whole population

 }
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survey}% at least one, from doc/KEYWORDS
