\name{svyglm}
\alias{svyglm}
\alias{svyglm.survey.design}
\alias{svyglm.svyrep.design}
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
Fit a generalised linear model to data from a complex survey design,
with inverse-probability weighting and design-based standard errors.
}
\usage{
\method{svyglm}{survey.design}(formula, design, subset=NULL, ...)
\method{svyglm}{svyrep.design}(formula, design, subset=NULL, ..., rho=NULL,
return.replicates=FALSE, na.action)
\method{summary}{svyglm}(object, correlation = FALSE, df.resid=NULL, ...) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula}
  \item{design}{Survey design from \code{\link{svydesign}} or \code{\link{svrepdesign}}. Must contain all variables
    in the formula}
  \item{subset}{Expression to select a subpopulation}
  \item{\dots}{Other arguments passed to \code{glm} or
    \code{summary.glm} }
  \item{rho}{For replicate BRR designs, to specify the parameter for
    Fay's variance method, giving weights of \code{rho} and \code{2-rho}}
  \item{return.replicates}{Return the replicates as a component of the result?}
  \item{object}{A \code{svyglm} object}
  \item{correlation}{Include the correlation matrix of parameters?}
  \item{na.action}{Handling of NAs}
  \item{df.resid}{Optional denominator degrees of freedom for Wald tests}
}
\details{
There is no \code{anova} method for \code{svyglm} as the models are not
fitted by maximum likelihood. The function \code{\link{regTermTest}} may
be useful for testing sets of regression terms.

If \code{df.resid} is not specified the df for the null model is
computed by \code{\link{degf}} and the residual df computed by
subtraction.  It's not that these are particularly good approximations
in a regression model but they are relatively standard. To get tests
based on a Normal distribution use \code{df.resid=Inf}.
}

 \value{ An object of class \code{svyglm}.  }

\author{Thomas Lumley}


\seealso{ \code{\link{svydesign}},
  \code{\link{svrepdesign}},\code{\link{as.svrepdesign}}, \code{\link{glm}},
  \code{\link{regTermTest}}}

\examples{

  data(api)

  glm(api00~ell+meals+mobility, data=apipop)

  dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  dclus2<-svydesign(id=~dnum+snum, weights=~pw, data=apiclus2)
  rstrat<-as.svrepdesign(dstrat)
  rclus2<-as.svrepdesign(dclus2)

  summary(svyglm(api00~ell+meals+mobility, design=dstrat))
  summary(svyglm(api00~ell+meals+mobility, design=dclus2))
  summary(svyglm(api00~ell+meals+mobility, design=rstrat))
  summary(svyglm(api00~ell+meals+mobility, design=rclus2))

  ## use quasibinomial, quasipoisson to avoid warning messages
  summary(svyglm(sch.wide~ell+meals+mobility, design=dstrat, family=quasibinomial()))

 }
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survey}% at least one, from doc/KEYWORDS
