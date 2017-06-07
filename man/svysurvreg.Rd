\name{svysurvreg}
\alias{svysurvreg}
\alias{svysurvreg.survey.design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit accelerated failure models to survey data
}
\description{
This function calls \code{survreg} from the 'survival' package to fit accelerated failure (accelerated life) models to complex survey data, and then computes correct standard errors by linearisation.  It has the same arguments as \code{survreg}, except that the second argument is \code{design} rather than \code{data}.  
}
\usage{
\method{svysurvreg}{survey.design}(formula, design, weights=NULL, subset=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Model formula
}
  \item{design}{
Survey design object, including two-phase designs
}
\item{weights}{
Additional weights to multiply by the sampling weights. No, I don't know why you'd want to do that.
}
  \item{subset}{
subset to use in fitting (if needed)
}

  \item{\dots}{
Other arguments of \code{survreg}
}
}
\value{
Object of class \code{svysurvreg}, with the same structure as a \code{survreg} object but with \code{NA} for the loglikelihood.
}

\note{
The \code{residuals} method is identical to that for \code{survreg} objects except the \code{weighted} option defaults to \code{TRUE}
}


\examples{

 data(pbc, package="survival")
 pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
 biasmodel<-glm(randomized~age*edema,data=pbc)
 pbc$randprob<-fitted(biasmodel)
 dpbc<-svydesign(id=~1, prob=~randprob, strata=~edema,
    data=subset(pbc,randomized))

 model <- svysurvreg(Surv(time, status>0)~bili+protime+albumin, design=dpbc, dist="weibull")
summary(model)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}% use one of  RShowDoc("KEYWORDS")
\keyword{survival}% __ONLY ONE__ keyword per line
