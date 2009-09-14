\name{regTermTest}
\alias{regTermTest}
\alias{print.regTermTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wald test for a term in a regression model}
\description{
  Provides Wald test and working likelihood ratio (Rao-Scott) test of the
  hypothesis that all coefficients associated with a particular
  regression term are zero (or have some other specified
  values). Particularly useful as a substitute for \code{\link{anova}}
  when not fitting by maximum likelihood.  The Wald tests use a chisquared or F distribution, the LRT uses a linear combination of chisquared distributions as in \code{\link{pchisqsum}}.
}
\usage{
regTermTest(model, test.terms, null=NULL,df=Inf, method=c("Wald","LRT"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A model object with \code{\link{coef}} and \code{\link{vcov}} methods}
  \item{test.terms}{Character string or one-sided formula giving name of
    term or terms to test}
  \item{null}{Null hypothesis values for parameters. Default is zeros}
  \item{df}{Denominator degrees of freedom for an F test. If
    \code{NULL} these are estimated from the model}
\item{method}{If \code{"Wald"}, the Wald-type test; if \code{"LRT"}
the Rao-Scott test based on the estimated log likelihood ratio}
}
\value{
  An object of class \code{regTermTest} or \code{regTermTestLRT}.
}
\references{
Rao, JNK, Scott, AJ (1984) "On Chi-squared Tests For Multiway Contingency Tables with Proportions Estimated From Survey Data" Annals of Statistics 12:46-60.
}
\seealso{\code{\link{anova}}, \code{\link{vcov}}, \code{\link{contrasts}},\code{\link{pchisqsum}}}
\examples{
 data(esoph)
 model1 <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp * 
     alcgp, data = esoph, family = binomial())
 anova(model1)

 regTermTest(model1,"tobgp")
 regTermTest(model1,"tobgp:alcgp")
 regTermTest(model1, ~alcgp+tobgp:alcgp)


 data(api)
 dclus2<-svydesign(id=~dnum+snum, weights=~pw, data=apiclus2)
 model2<-svyglm(I(sch.wide=="Yes")~ell+meals+mobility, design=dclus2, family=quasibinomial())
 regTermTest(model2, ~ell)
 regTermTest(model2, ~ell,df=NULL)
 regTermTest(model2, ~ell, method="LRT")
 regTermTest(model2, ~ell+meals, method="LRT")
}
\keyword{regression}% at least one, from doc/KEYWORDS

