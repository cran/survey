\name{svyglm}
\alias{svyglm}
\alias{svrepglm}
\alias{print.svyglm}
\alias{summary.svyglm}
\alias{summary.svrepglm}
\alias{print.summary.svyglm}
\alias{vcov.svyglm}
\alias{residuals.svyglm}
\alias{residuals.svrepglm}
\alias{coef.svyglm}
\alias{extractAIC.svyglm}
\alias{extractAIC.svrepglm}
\alias{logLik.svyglm}
\alias{logLik.svrepglm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Survey-weighted generalised linear models.}
\description{
Fit a generalised linear model to data from a complex survey design, with inverse-probability weighting and with standard errors corrected for cluster sampling.
}
\usage{
svyglm(formula, design, subset=NULL, ...)
svrepglm(formula, design, subset=NULL, ..., rho=NULL,
return.replicates=FALSE, na.action)
\method{summary}{svyglm}(object, correlation = FALSE,  ...) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula}
  \item{design}{Survey design from \code{\link{svydesign}} or \code{\link{svrepdesign}}. Must contain all variables
    in the formula}
  \item{subset}{Expression to select a subpopulation}
  \item{\dots}{Other arguments passed to \code{glm} or
    \code{summary.glm} }
  \item{rho}{For replicate BRR designs, to specify the paramter for
    Fay's variance method}
  \item{return.replicates}{Return the replicates as a component of the result?}
  \item{object}{A \code{svyglm} object}
  \item{correlation}{Include the correlation matrix of parameters?}
  \item{na.action}{Handling of NAs}
}
\details{
In \code{svyglm}, standard errors for cluster-sampled designs are computed using a
linearisation estimate (in the absence of strata this is equivalent to
the Huber/White sandwich formula in GEEs).  Most of these computations
are done in \code{\link{svyCprod}}.  In \code{svrepglm}, replicate
weight methods are used.

There is no \code{anova} method for \code{svyglm} as the models are not
fitted by maximum likelihood. The function \code{\link{regTermTest}} may
be useful for testing sets of regression terms.
}
\value{
  An object of class \code{svyglm}.
}

\author{Thomas Lumley}


\seealso{ \code{\link{svydesign}},
  \code{\link{svrepdesign}},\code{\link{as.svrepdesign}}, \code{\link{glm}},
  \code{\link{svyCprod}}, \code{\link{svy.varcoef}},\code{\link{regTermTest}}}

\examples{
  ## Independent sampling
  df<-data.frame(x=rnorm(1000),z=rep(0:4,200))
  df$y<-with(df, 3+3*x*z)
  df$p<-with(df, exp(x)/(1+exp(x)))
  mpop<-lm(y~x+z, data=df) 
  xi<-rbinom(1000,1,df$p)
  sdf<-df[xi==1,]

  dxi<-svydesign(~0,~p,data=sdf)
  m<-svyglm(y~x+z,family="gaussian",design=dxi)
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
  drep<-as.svrepdesign(dxi)
  m<-svyglm(y~x+z,family="gaussian",design=dxi)
  mr<-svrepglm(y~x+z, family="gaussian", design=drep)
  m1<-lm(y~x+z,data=sdf)
 
  summary(m1) ##wrong
  summary(m)  ##right
  summary(mpop) ##whole population

  ## subsets
  msub<-svyglm(y~x+z,family="gaussian",design=dxi,subset=x>1)
  summary(msub)
  subdxi<-subset(dxi,x>1)
  msub<-svyglm(y~x+z,family="gaussian",design=subdxi)
  summary(msub)

 }
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survey}% at least one, from doc/KEYWORDS
