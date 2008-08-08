\name{svyolr}
\alias{svyolr}
\alias{svyolr.survey.design2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Proportional odds and related models }
\description{
  Fits cumulative link models: proportional odds, probit, complementary
  log-log, and cauchit. 
}
\usage{
svyolr(formula, design, ...)
\method{svyolr}{survey.design2}(formula, design, start, ..., na.action = na.omit, method = c("logistic", 
    "probit", "cloglog", "cauchit")) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Formula: the response must be a factor with at least
    three levels}
  \item{design}{survey design object }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
  \item{start}{Optional starting values for optimization}
  \item{na.action}{handling of missing values}
  \item{method}{Link function}
}
\value{
  An object of class \code{svyolr}
}

\author{The code is based closely on polr() from the MASS package of
  Venables and Ripley.}

\seealso{\code{\link{svyglm}} }
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
dclus1<-update(dclus1, mealcat=cut(meals,c(0,25,50,75,100)))

svyolr(mealcat~avg.ed+mobility+stype, design=dclus1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

