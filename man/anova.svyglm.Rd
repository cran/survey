\name{anova.svyglm}
\alias{anova.svyglm}
\alias{AIC.svyglm}
\alias{BIC.svyglm}
\alias{extractAIC.svyglm}
\alias{extractAIC.svrepglm}
\alias{anova.svycoxph}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model comparison for glms.
}
\description{
A method for the \code{\link{anova}} function, for use on
\code{\link{svyglm}} and \code{\link{svycoxph}} objects.  With a single model argument it produces a sequential anova table, with two arguments it compares the two models.
}
\usage{
\method{anova}{svyglm}(object, object2 = NULL, test = c("F", "Chisq"), 
 method = c("LRT", "Wald"), tolerance = 1e-05, ..., force = FALSE)
\method{anova}{svycoxph}(object, object2=NULL,test=c("F","Chisq"),
 method=c("LRT","Wald"),tolerance=1e-5,...,force=FALSE)
\method{AIC}{svyglm}(object,...,k=2, null_has_intercept=TRUE)
\method{BIC}{svyglm}(object,...,maximal)
\method{extractAIC}{svyglm}(fit,scale,k=2,..., null_has_intercept=TRUE)
\method{extractAIC}{svrepglm}(fit,scale,k=2,..., null_has_intercept=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object,fit}{
  A \code{\link{svyglm}} or\code{\link{svycoxph}}  object.
}
  \item{object2}{
  Optionally, another \code{\link{svyglm}} or \code{\link{svycoxph}} object.
}
  \item{test}{
  Use (linear combination of) F or chi-squared distributions for p-values. F is usually preferable.
}
  \item{method}{
  Use weighted deviance difference (LRT) or Wald tests to compare models
}
  \item{tolerance}{
  For models that are not symbolically nested, the tolerance for deciding that a term is common to the models.
}
  \item{\dots}{
  For \code{AIC} and \code{BIC}, optionally more \code{svyglm} objects
}
\item{scale}{not used}
\item{null_has_intercept}{Does the null model for AIC have an intercept or
  not? Must be \code{FALSE} if any of the models are intercept-only. }
  \item{force}{
   Force the tests to be done by explicit projection even if the models
   are symbolically nested (eg, for debugging)
}
\item{maximal}{A \code{svyglm} model that \code{object} (and \dots if supplied) are nested in.}
\item{k}{Multiplier for effective df in AIC. Usually 2. There is no choice of \code{k} that will give BIC}
}
\details{
The reference distribution for the LRT depends on the misspecification effects for the parameters being tested (Rao and Scott, 1984). If the models are symbolically nested, so that the relevant parameters can be identified just by manipulating the model formulas, \code{anova} is equivalent to \code{\link{regTermTest}}.If the models are nested but not symbolically nested, more computation using the design matrices is needed to determine the projection matrix on to the parameters being tested.   In the examples below, \code{model1} and \code{model2} are symbolically nested in \code{model0} because \code{model0} can be obtained just by deleting terms from the formulas.  On the other hand, \code{model2} is nested in \code{model1} but not symbolically nested: knowing that the model is nested requires knowing what design matrix columns are produced by \code{stype} and \code{as.numeric(stype)}.  Other typical examples of models that are nested but not symbolically nested are linear and spline models for a continuous covariate, or models with categorical versions of a variable at different resolutions (eg, smoking yes/no or smoking never/former/current). 

A saddlepoint approximation is used for the LRT with numerator df greater than 1.

\code{AIC} is defined using the Rao-Scott approximation to the weighted
loglikelihood (Lumley and Scott, 2015). It replaces the usual penalty term p, which is the null expectation of the log likelihood ratio, by the trace of the generalised design effect matrix, which is the expectation under complex sampling. For computational reasons everything is scaled so the weights sum to the sample size. 

\code{BIC} is a BIC for the (approximate) multivariate Gaussian models
on regression coefficients from the maximal model implied by each
submodel (ie, the models that say some coefficients in the maximal model
are zero) (Lumley and Scott, 2015). It corresponds to comparing the models with a Wald test and replacing the sample size in the penalty by an effective sample size.
For computational reasons, the models must not only be nested, the names of the coefficients must match.

\code{extractAIC} for a model with a Gaussian link uses the actual AIC based on maximum likelihood estimation of the variance parameter as well as the regression parameters.

}
\value{
Object of class \code{seqanova.svyglm} if one model is given, otherwise of class \code{regTermTest} or \code{regTermTestLRT}
}
\note{
At the moment, \code{AIC} works only for models including an intercept.
  }
\references{

Rao, JNK, Scott, AJ (1984) "On Chi-squared Tests For Multiway Contingency Tables with Proportions Estimated From Survey Data" Annals of Statistics 12:46-60.

Lumley, T., & Scott, A. (2014). "Tests for Regression Models Fitted to Survey Data". Australian and New Zealand Journal of Statistics, 56 (1), 1-14. 

Lumley T, Scott AJ (2015)
"AIC and BIC for modelling with complex survey data" J Surv Stat
Methodol 3 (1): 1-18. 

}

\seealso{
\code{\link{regTermTest}}, \code{\link{pchisqsum}}
}
\examples{
data(api)
dclus2<-svydesign(id=~dnum+snum, weights=~pw, data=apiclus2)

model0<-svyglm(I(sch.wide=="Yes")~ell+meals+mobility, design=dclus2, family=quasibinomial())
model1<-svyglm(I(sch.wide=="Yes")~ell+meals+mobility+as.numeric(stype), 
     design=dclus2, family=quasibinomial())
model2<-svyglm(I(sch.wide=="Yes")~ell+meals+mobility+stype, design=dclus2, family=quasibinomial())

anova(model2)	
anova(model0,model2)					     		    
anova(model1, model2)

anova(model1, model2, method="Wald")

AIC(model0,model1, model2)
BIC(model0, model2,maximal=model2)

## AIC for linear model is different because it considers the variance
## parameter

model0<-svyglm(api00~ell+meals+mobility, design=dclus2)
model1<-svyglm(api00~ell+meals+mobility+as.numeric(stype), 
     design=dclus2)
model2<-svyglm(api00~ell+meals+mobility+stype, design=dclus2)
modelnull<-svyglm(api00~1, design=dclus2)

AIC(model0, model1, model2)

AIC(model0, model1, model2,modelnull, null_has_intercept=FALSE)

## from ?twophase
data(nwtco)
dcchs<-twophase(id=list(~seqno,~seqno), strata=list(NULL,~rel),
        subset=~I(in.subcohort | rel), data=nwtco)
a<-svycoxph(Surv(edrel,rel)~factor(stage)+factor(histol)+I(age/12), design=dcchs)
b<-update(a, .~.-I(age/12)+poly(age,3))
## not symbolically nested models
anova(a,b)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{regression}% __ONLY ONE__ keyword per line
