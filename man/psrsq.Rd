\name{psrsq}
\alias{psrsq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pseudo-Rsquareds
}
\description{
Compute the Nagelkerke and Cox--Snell pseudo-rsquared statistics, primarily for logistic regression. A generic function with methods for \code{glm} and \code{\link{svyglm}}.  The method for \code{svyglm} objects uses the design-based estimators described by Lumley (2017)
}
\usage{
psrsq(object, method = c("Cox-Snell", "Nagelkerke"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A regression model (\code{glm} or \code{svyglm})
}
  \item{method}{
Which statistic to compute
}
  \item{\dots}{
For future expansion
}
}

\value{
Numeric value
}
\references{
Lumley T (2017) "Pseudo-R2 statistics under complex sampling" Australian and New Zealand Journal of Statistics DOI: 10.1111/anzs.12187 (preprint: \url{https://arxiv.org/abs/1701.07745})
}

\seealso{
\code{\link{AIC.svyglm}}
}
\examples{
data(api)
dclus2<-svydesign(id=~dnum+snum, weights=~pw, data=apiclus2)

model1<-svyglm(I(sch.wide=="Yes")~ell+meals+mobility+as.numeric(stype), 
     design=dclus2, family=quasibinomial())

psrsq(model1, type="Nagelkerke")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey }% use one of  RShowDoc("KEYWORDS")
\keyword{regression }% __ONLY ONE__ keyword per line
