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
\method{confint}{svyglm}(object, parm, level = 0.95, method = c("Wald", "likelihood"), ddf = Inf, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{svyglm} object}
  \item{parm}{numeric or character vector indicating which parameters to
  construct intervals for.}
  \item{level}{desired coverage}
  \item{method}{See description above }
  \item{ddf}{Denominator degrees of freedom for \code{"likelihood"}
    method, to use a t distribution rather than norma. If \code{NULL},
    use \code{object$df.residual}}
  \item{\dots}{for future expansion}
}

\value{
  A matrix of confidence intervals
}
\references{
J. N. K. Rao and Alistair J. Scott (1984) On Chi-squared Tests For
  Multiway Contigency Tables with Proportions Estimated From Survey
  Data. Annals of Statistics 12:46-60
}

\seealso{\code{\link{confint}} }
\examples{
data(api)
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)

m<-svyglm(I(comp.imp=="Yes")~stype*emer+ell, design=dclus2, family=quasibinomial)
confint(m)
confint(m, method="like",ddf=NULL, parm=c("ell","emer"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

