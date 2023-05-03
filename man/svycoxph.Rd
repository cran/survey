\name{svycoxph}
\alias{svycoxph}
\alias{svycoxph.survey.design2}
\alias{svycoxph.survey.design}
\alias{svycoxph.svyrep.design}
\alias{predict.svycoxph}
%\alias{print.svycoxph}
%\alias{model.frame.svycoxph}
%\alias{summary.svycoxph}
%\alias{anova.svycoxph}
\alias{extractAIC.svycoxph}
\alias{AIC.svycoxph}
%\alias{survfit.svycoxph}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Survey-weighted Cox models.}
\description{
Fit a proportional hazards model to data from a complex survey design.
}
\usage{
svycoxph(formula, design,subset=NULL, rescale=TRUE, ...)
\method{predict}{svycoxph}(object, newdata, se=FALSE,
    type=c("lp", "risk", "terms","curve"),...)
\method{AIC}{svycoxph}(object, ..., k = 2)    
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula. Any \code{cluster()} terms will be ignored.}
  \item{design}{ \code{survey.design} object. Must contain all variables
    in the formula}
  \item{subset}{Expression to select a subpopulation}
  \item{rescale}{Rescale weights to improve numerical stability}
  \item{object}{A \code{svycoxph} object}
  \item{newdata}{New data for prediction}
  \item{se}{Compute standard errors? This takes a lot of memory for
    \code{type="curve"}}
  \item{type}{"curve" does predicted survival curves. The other values
    are passed to \code{predict.coxph()}}
  \item{\dots}{For \code{AIC}, more models to compare the AIC of. For \code{svycoxph},
  other arguments passed to \code{coxph}. }
  \item{k}{The penalty per parameter that would be used under independent sampling: AIC has \code{k=2}}
}
\details{
The main difference between \code{svycoxph} function and the \code{robust=TRUE}
option to \code{\link{coxph}} in the
survival package is that this function accounts for the reduction in
variance from stratified sampling and the increase in variance from
having only a small number of clusters.

Note that \code{strata} terms in the model formula describe subsets that
have a separate baseline hazard function and need not have anything to
do with the stratification of the sampling.

The \code{AIC} method uses the same approach as \code{\link{AIC.svyglm}},
though the relevance of the criterion this optimises is a bit less clear
than for generalised linear models. 

The standard errors for predicted survival curves are available only by linearization, not
  by replicate weights (at the moment). Use
  \code{\link{withReplicates}} to get standard errors with replicate
  weights. Predicted survival curves are not available for stratified
  Cox models.

The standard errors use the delta-method approach of Williams (1995)
  for the Nelson-Aalen estimator, modified to handle the Cox model
  following Tsiatis (1981). The standard errors agree closely with
  \code{survfit.coxph} for independent sampling when the model fits
  well, but are larger when the model fits poorly.  I believe the
  standard errors are equivalent to those of Lin (2000), but I don't
  know of any implementation that would allow a check.
  }
\value{
  An object of class \code{svycoxph} for \code{svycoxph}, an object of
  class \code{svykm} or \code{svykmlist} for \code{predict(,type="curve")}.
}
\section{Warning}{
The standard error calculation for survival curves uses memory
  proportional to the sample size times the square of the number of events.
  }
\author{Thomas Lumley}

\references{Binder DA. (1992) Fitting Cox's proportional hazards models
  from survey data. Biometrika 79: 139-147

Lin D-Y (2000) On fitting Cox's proportional hazards model to survey data. Biometrika 87: 37-47

Tsiatis AA (1981) A Large Sample Study of Cox's Regression Model. Annals
of Statistics 9(1) 93-108

Williams RL (1995) "Product-Limit Survival Functions with Correlated 
Survival Times" Lifetime Data Analysis 1: 171--186

}
\seealso{  \code{\link{coxph}}, \code{\link{predict.coxph}}

  \code{\link{svykm}} for estimation of Kaplan-Meier survival curves and
  for methods that operate on survival curves.

\code{\link{regTermTest}} for Wald and (Rao-Scott) likelihood ratio tests for one or more parameters.
}

\examples{
## Somewhat unrealistic example of nonresponse bias.
data(pbc, package="survival")

pbc$randomized<-with(pbc, !is.na(trt) & trt>0)
biasmodel<-glm(randomized~age*edema,data=pbc,family=binomial)
pbc$randprob<-fitted(biasmodel)
if (is.null(pbc$albumin)) pbc$albumin<-pbc$alb ##pre2.9.0

dpbc<-svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,randomized))
rpbc<-as.svrepdesign(dpbc)

(model<-svycoxph(Surv(time,status>0)~log(bili)+protime+albumin,design=dpbc))

svycoxph(Surv(time,status>0)~log(bili)+protime+albumin,design=rpbc)

s<-predict(model,se=TRUE, type="curve",
     newdata=data.frame(bili=c(3,9), protime=c(10,10), albumin=c(3.5,3.5)))
plot(s[[1]],ci=TRUE,col="sienna")
lines(s[[2]], ci=TRUE,col="royalblue")
quantile(s[[1]], ci=TRUE)
confint(s[[2]], parm=365*(1:5))
}
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survival}% at least one, from doc/KEYWORDS
\keyword{survey}% at least one, from doc/KEYWORDS
