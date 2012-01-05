\name{anova.svyglm}
\alias{anova.svyglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model comparison for glms.
}
\description{
A method for the \code{\link{anova}} function, for use on \code{\link{svyglm}} objects.  With a single model argument it produces a sequential anova table, with two arguments it compares the two models.
}
\usage{
\method{anova}{svyglm}(object, object2 = NULL, test = c("F", "Chisq"), method = c("LRT", "Wald"), tolerance = 1e-05, ..., force = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A \code{\link{svyglm}} object.
}
  \item{object2}{
  Optionally, another \code{\link{svyglm}} object.
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
  For future expansion
}
  \item{force}{
   Force the tests to be done by explicit projection even if the models are symbolically nested (for debugging)
}
}
\details{
The reference distribution for the LRT depends on the misspecification effects for the parameters being tested (Rao and Scott, 1984). If the models are symbolically nested, so that the relevant parameters can be identified just by manipulating the model formulas, \code{anova} is equivalent to \code{\link{regTermTest}}. If the models are nested but not symbolically nested, more computation using the design matrices is needed to determine the projection matrix on to the parameters being tested.  Typical examples of models that are nested but not symbolically nested are linear and spline models for a continuous covariate or linear and saturated models for a factor. 

The saddlepoint approximation is used for the LRT with numerator df greater than 1.
}
\value{
Object of class \code{seqanova.svyglm} if one model is given, otherwise of class \code{regTermTest} or \code{regTermTestLRT}
}
\references{

Rao, JNK, Scott, AJ (1984) "On Chi-squared Tests For Multiway Contingency Tables with Proportions Estimated From Survey Data" Annals of Statistics 12:46-60.
}

\seealso{
\code{\link{regTermTest}}, \code{\link{pchisqsum}}
}
\examples{
data(api)
 dclus2<-svydesign(id=~dnum+snum, weights=~pw, data=apiclus2)

 model1<-svyglm(I(sch.wide=="Yes")~ell+meals+mobility+as.numeric(stype), design=dclus2, family=quasibinomial())
 model2<-svyglm(I(sch.wide=="Yes")~ell+meals+mobility+stype, design=dclus2, family=quasibinomial())

 anova(model2)						     		    
 anova(model1, model2)

 anova(model1, model2, method="Wald")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{regression}% __ONLY ONE__ keyword per line
