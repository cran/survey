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

  data(api)

  glm(api00~ell+meals+mobility, data=apipop)

  dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  dclus2<-svydesign(id=~dnum+snum, weights=~pw, data=apiclus2)
  rstrat<-as.svrepdesign(dstrat)
  rclus2<-as.svrepdesign(dclus2)

  summary(svyglm(api00~ell+meals+mobility, design=dstrat))
  summary(svyglm(api00~ell+meals+mobility, design=dclus2))
  summary(svrepglm(api00~ell+meals+mobility, design=rstrat))
  summary(svrepglm(api00~ell+meals+mobility, design=rclus2))

  ## use quasibinomial, quasipoisson to avoid warning messages
  summary(svyglm(sch.wide~ell+meals+mobility, design=dstrat, family=quasibinomial()))

 }
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survey}% at least one, from doc/KEYWORDS
