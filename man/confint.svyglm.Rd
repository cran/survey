\name{confint.svyglm}
\alias{confint.svyglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence intervals for regression parameters }
\description{
Computes confidence intervals for regression parameters in
\code{\link{svyglm}} objects. The default is a Wald-type confidence
interval, adding and subtracting a multiple of the standard error. The
\code{method="likelihood"} is an interval based on inverting the Rao-Scott
likelihood ratio test. That is, it is an interval where the working
model deviance is lower than the threshold for the Rao-Scott test at the
specified level.
}
\usage{
\method{confint}{svyglm}(object, parm, level = 0.95, method = c("Wald", "likelihood"), ddf = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{\link{svyglm}} object}
  \item{parm}{numeric or character vector indicating which parameters to
  construct intervals for.}
  \item{level}{desired coverage}
  \item{method}{See description above }
  \item{ddf}{Denominator degrees of freedom for \code{"likelihood"}
    method, to use a t distribution rather than normal. If \code{NULL},
    use \code{object$df.residual}
    for Taylor-series standard errors, or \code{object$df.coef}
    for Bell-McCaffrey standard errors with adjusted degrees of freedom.}
  \item{\dots}{for future expansion}
}

\value{
  A matrix of confidence intervals, 
  possibly with additional attributes \code{levels} and/or \code{degf}.
}
\references{
J. N. K. Rao and Alastair J. Scott (1984). On Chi-squared Tests For
  Multiway Contigency Tables with Proportions Estimated From Survey
  Data. Annals of Statistics 12:46-60.
  

Robert M. Bell and Daniel F. McCaffrey (2002).
  Bias Reduction in Standard Errors for Linear Regression with
  Multi-Stage Samples. Survey Methodology 28 (2), 169-181.
  \url{https://www150.statcan.gc.ca/n1/pub/12-001-x/2002002/article/9058-eng.pdf}
}

\seealso{\code{\link{confint}} }
\examples{
data(api)
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)

m<-svyglm(I(comp.imp=="Yes")~stype*emer+ell, design=dclus2, family=quasibinomial)
confint(m)
confint(m, method="like",ddf=NULL, parm=c("ell","emer"))

m2<-svyglm(I(comp.imp=="Yes")~stype*emer+ell, design=dclus2, family=quasibinomial,
    std.errors="Bell-McCaffrey-2", degf=TRUE)
confint(m2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
