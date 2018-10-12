\name{svyglm}
\alias{svyglm}
\alias{svyglm.survey.design}
\alias{svyglm.svyrep.design}
\alias{summary.svyglm}
\alias{summary.svrepglm}
\alias{vcov.svyglm}
\alias{residuals.svyglm}
\alias{residuals.svrepglm}
\alias{predict.svyglm}
\alias{predict.svrepglm}
\alias{coef.svyglm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Survey-weighted generalised linear models.}
\description{
Fit a generalised linear model to data from a complex survey design,
with inverse-probability weighting and design-based standard errors.
}
\usage{
\method{svyglm}{survey.design}(formula, design, subset=NULL,
family=stats::gaussian(),start=NULL, rescale=TRUE, ...)
\method{svyglm}{svyrep.design}(formula, design, subset=NULL,
family=stats::gaussian(),start=NULL, rescale=NULL, ..., rho=NULL,
return.replicates=FALSE, na.action,multicore=getOption("survey.multicore"))
\method{summary}{svyglm}(object, correlation = FALSE, df.resid=NULL,
...)
\method{predict}{svyglm}(object,newdata=NULL,total=NULL,
                         type=c("link","response","terms"),
                         se.fit=(type != "terms"),vcov=FALSE,...)
\method{predict}{svrepglm}(object,newdata=NULL,total=NULL,
                         type=c("link","response","terms"),
                         se.fit=(type != "terms"),vcov=FALSE,
                         return.replicates=!is.null(object$replicates),...)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula}
  \item{design}{Survey design from \code{\link{svydesign}} or \code{\link{svrepdesign}}. Must contain all variables
    in the formula}
  \item{subset}{Expression to select a subpopulation}
  \item{family}{\code{family} object for \code{glm}}
  \item{start}{Starting values for the coefficients (needed for some
    uncommon link/family combinations)}
  \item{rescale}{Rescaling of weights, to improve numerical stability. The default
    rescales weights to sum to the sample size. Use \code{FALSE} to not
    rescale weights. For replicate-weight designs, use \code{TRUE} to
    rescale weights to sum to 1, as was the case before version 3.34.}
  \item{\dots}{Other arguments passed to \code{glm} or
    \code{summary.glm} }
  \item{rho}{For replicate BRR designs, to specify the parameter for
    Fay's variance method, giving weights of \code{rho} and \code{2-rho}}
  \item{return.replicates}{Return the replicates as the \code{replicates} component of the
    result? (for \code{predict}, only possible if they
    were computed in the \code{svyglm} fit)}
  \item{object}{A \code{svyglm} object}
  \item{correlation}{Include the correlation matrix of parameters?}
  \item{na.action}{Handling of NAs}
  \item{multicore}{Use the \code{multicore} package to distribute
    replicates across processors?}	      
  \item{df.resid}{Optional denominator degrees of freedom for Wald
    tests}
  \item{newdata}{new data frame for prediction}
  \item{total}{population size when predicting population total}
  \item{type}{linear predictor (\code{link}) or response}
  \item{se.fit}{if \code{TRUE}, return variances of predictions}
  \item{vcov}{if \code{TRUE} and \code{se=TRUE} return full
    variance-covariance matrix of predictions}
}
\details{
For binomial and Poisson families use \code{family=quasibinomial()}
and \code{family=quasipoisson()} to avoid a warning about non-integer
numbers of successes.  The `quasi' versions of the family objects give
the same point estimates and standard errors and do not give the
warning.

If \code{df.resid} is not specified the df for the null model is
computed by \code{\link{degf}} and the residual df computed by
subtraction.  This is recommended by Korn and Graubard and is correct
for PSU-level covariates but is potentially very conservative for
individual-level covariates. To get tests based on a Normal distribution
use \code{df.resid=Inf}, and to use number of PSUs-number of strata,
specify \code{df.resid=degf(design)}.

Parallel processing with \code{multicore=TRUE} is helpful only for
fairly large data sets and on computers with sufficient memory. It may
be incompatible with GUIs, although the Mac Aqua GUI appears to be safe.

\code{predict} gives fitted values and sampling variability for specific new
values of covariates.  When \code{newdata} are the population mean it
gives the regression estimator of the mean, and when \code{newdata} are
the population totals and \code{total} is specified it gives the
regression estimator of the population total.  Regression estimators of
mean and total can also be obtained with \code{\link{calibrate}}.


}

 \value{ \code{svyglm} returns an object of class \code{svyglm}.  The
 \code{predict} method returns an object of class \code{svystat}}

\author{Thomas Lumley}


\seealso{ 
  \code{\link{glm}}, which is used to do most of the work.
  
  \code{\link{regTermTest}}, for multiparameter tests

  \code{\link{calibrate}}, for an alternative way to specify regression
  estimators of population totals or means

  \code{\link{svyttest}} for one-sample and two-sample t-tests.
}

\references{

Lumley T, Scott A (2017) "Fitting Regression Models to Survey Data"
Statistical Science 32: 265-278

  }

\examples{

  data(api)


  dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  dclus2<-svydesign(id=~dnum+snum, weights=~pw, data=apiclus2)
  rstrat<-as.svrepdesign(dstrat)
  rclus2<-as.svrepdesign(dclus2)

  summary(svyglm(api00~ell+meals+mobility, design=dstrat))
  summary(svyglm(api00~ell+meals+mobility, design=dclus2))
  summary(svyglm(api00~ell+meals+mobility, design=rstrat))
  summary(svyglm(api00~ell+meals+mobility, design=rclus2))

  ## use quasibinomial, quasipoisson to avoid warning messages
  summary(svyglm(sch.wide~ell+meals+mobility, design=dstrat,
        family=quasibinomial()))


  ## Compare regression and ratio estimation of totals
  api.ratio <- svyratio(~api.stu,~enroll, design=dstrat)
  pop<-data.frame(enroll=sum(apipop$enroll, na.rm=TRUE))
  npop <- nrow(apipop)
  predict(api.ratio, pop$enroll)

  ## regression estimator is less efficient
  api.reg <- svyglm(api.stu~enroll, design=dstrat)
  predict(api.reg, newdata=pop, total=npop)
  ## same as calibration estimator
  svytotal(~api.stu, calibrate(dstrat, ~enroll, pop=c(npop, pop$enroll)))

  ## svyglm can also reproduce the ratio estimator
  api.reg2 <- svyglm(api.stu~enroll-1, design=dstrat,
                    family=quasi(link="identity",var="mu"))
  predict(api.reg2, newdata=pop, total=npop)

  ## higher efficiency by modelling variance better
  api.reg3 <- svyglm(api.stu~enroll-1, design=dstrat,
                    family=quasi(link="identity",var="mu^3"))
  predict(api.reg3, newdata=pop, total=npop)
  ## true value
  sum(apipop$api.stu)

 }
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survey}% at least one, from doc/KEYWORDS
