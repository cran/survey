\name{svycoxph}
\alias{svycoxph}
\alias{print.svycoxph}
\alias{model.frame.svycoxph}
\alias{summary.svycoxph}
\alias{anova.svycoxph}
\alias{extractAIC.svycoxph}
\alias{survfit.svycoxph}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Survey-weighted generalised linear models.}
\description{
Fit a generalised linear model to data from a complex survey design, with inverse-probability weighting and with standard errors corrected for cluster sampling.
}
\usage{
svycoxph(formula, design,subset=NULL, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula. Any \code{cluster()} terms will be ignored.}
  \item{design}{ \code{survey.design} object. Must contain all variables
    in the formula}
  \item{subset}{Expression to select a subpopulation}
  \item{\dots}{Other arguments passed to \code{coxph}. }
}
\details{
The main difference between this function and \code{\link{coxph}} in the
survival package is that this function accounts for the reduction in
variance from stratified sampling and the increase in variance from
having only a small number of clusters.

Note that \code{strata} terms in the model formula describe subsets that
have a separate baseline hazard function and need not have anything to
do with the stratification of the sampling.
  }
\value{
  An object of class \code{svycoxph}.
}

\author{Thomas Lumley}


\seealso{ \code{\link{svydesign}}, \code{\link{coxph}},
  \code{\link{svyCprod}}}

\examples{

 }
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{survival}% at least one, from doc/KEYWORDS
\keyword{survey}% at least one, from doc/KEYWORDS
