\name{svyolr}
\alias{svyolr}
\alias{svyolr.survey.design2}
\alias{svyolr.svyrep.design}
\alias{predict.svyolr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Proportional odds and related models }
\description{
  Fits cumulative link models: proportional odds, probit, complementary
  log-log, and cauchit.
}
\usage{
svyolr(formula, design, ...)
\method{svyolr}{survey.design2}(formula, design, start, subset=NULL,...,
    na.action = na.omit,method = c("logistic", "probit", "cloglog", "cauchit"))
\method{svyolr}{svyrep.design}(formula,design,subset=NULL,...,return.replicates=FALSE, 
    multicore=getOption("survey.multicore"))
\method{predict}{svyolr}(object, newdata, type = c("class", "probs"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Formula: the response must be a factor with at least
    three levels}
  \item{design}{survey design object }
  \item{subset}{subset of the design to use; \code{NULL} for all of it}
  \item{\dots}{dots}
  \item{start}{Optional starting values for optimization}
  \item{na.action}{handling of missing values}
  \item{multicore}{Use \code{multicore} package to distribute computation of replicates across multiple
    processors?}			   
  \item{method}{Link function}
  \item{return.replicates}{return the individual replicate-weight
    estimates}
  \item{object}{object of class \code{svyolr}}
  \item{newdata}{new data for predictions}
  \item{type}{return vector of most likely class or matrix of probabilities}
}
\value{
  An object of class \code{svyolr}
}

\author{The code is based closely on polr() from the MASS package of
  Venables and Ripley.}

\seealso{\code{\link{svyglm}}, \code{\link{regTermTest}} }
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
dclus1<-update(dclus1, mealcat=cut(meals,c(0,25,50,75,100)))

m<-svyolr(mealcat~avg.ed+mobility+stype, design=dclus1)
m

## Use regTermTest for testing multiple parameters
regTermTest(m, ~avg.ed+stype, method="LRT")

## predictions
summary(predict(m, newdata=apiclus2))
summary(predict(m, newdata=apiclus2, type="probs"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

