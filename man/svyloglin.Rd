\name{svyloglin}
\alias{svyloglin}
\alias{anova.svyloglin}
\alias{update.svyloglin}
\alias{coef.svyloglin}
\alias{print.anova.svyloglin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Loglinear models }
\description{
Fit and compare hierarchical loglinear models for complex survey data.
}
\usage{
svyloglin(formula, design, ...)
\method{update}{svyloglin}(object,formula,...)
\method{anova}{svyloglin}(object,object1,...,integrate=FALSE)
\method{print}{anova.svyloglin}(x,pval=c("F","saddlepoint","lincom","chisq"),...)
\method{coef}{svyloglin}(object,...,intercept=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula}
  \item{design}{survey design object}
  \item{object,object1}{loglinear model from \code{svyloglin}}
  \item{pval}{p-value approximation: see Details}
  \item{integrate}{Compute the exact asymptotic p-value (slow)?}
  \item{\dots}{not used }
  \item{intercept}{Report the intercept?}
  \item{x}{anova object}
}
\details{
The loglinear model is fitted to a multiway table with probabilities
estimated by \code{\link{svymean}} and with the sample size equal to the
observed sample size, treating the resulting table as if it came from iid
multinomial sampling, as described by Rao and Scott. The
variance-covariance matrix does not include the intercept term, and so
by default neither does the \code{coef} method.  A Newton-Raphson
algorithm is used, rather than iterative proportional fitting, so
starting values are not needed.

The \code{anova} method computes the quantities that would be the score
(Pearson) and likelihood ratio chi-squared statistics if the data were
an iid sample. It computes four p-values for each of these, based on the
exact asymptotic distribution (see \code{\link{pchisqsum}}), a
saddlepoint approximateion to this distribution, a scaled
chi-squared distribution, and a scaled F-distribution.  When testing the
two-way interaction model against the main-effects model in a two-way
table the score statistic and p-values match the Rao-Scott tests
computed by \code{\link{svychisq}}.

The \code{anova} method can only compare two models if they are for
exactly the same multiway table (same variables and same order). The
\code{update} method will help with this. It is also much faster to use
\code{update} than \code{svyloglin} for a large data set: its time
complexity depends only on the size of the model, not on the size of the
data set.

 It is not possible to fit a model using a variable created inline, eg
\code{I(x<10)}, since the multiway table is based on all variables used
in the formula. 
}
\value{
  Object of class \code{"svyloglin"}
}
\references{
Rao, JNK, Scott, AJ (1984) "On Chi-squared Tests For Multiway Contingency Tables with Proportions Estimated From Survey Data" Annals of Statistics 12:46-60.
}

\seealso{\code{\link{svychisq}}, \code{\link{svyglm}},\code{\link{pchisqsum}}}
\examples{
 data(api)
 dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
 a<-svyloglin(~stype+comp.imp,dclus1)
 b<-update(a,~.^2)
 an<-anova(a,b)
 an
 print(an, pval="saddlepoint")

 ## Wald test
 regTermTest(b, ~stype:comp.imp)

 ## linear-by-linear association
 d<-update(a,~.+as.numeric(stype):as.numeric(comp.imp))
 an1<-anova(a,d)
 an1


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

