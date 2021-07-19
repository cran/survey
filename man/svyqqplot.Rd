\name{svyqqplot}
\alias{svyqqplot}
\alias{svyqqmath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile-quantile plots for survey data
}
\description{
Quantile-quantile plots either against a specified distribution function or comparing two variables from the same or different designs.
}
\usage{
svyqqplot(formula, design, designx = NULL, na.rm = TRUE, qrule = "hf8",
    xlab = NULL, ylab = NULL, ...)
svyqqmath(x, design, null=qnorm, na.rm=TRUE, xlab="Expected",ylab="Observed",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,formula}{
A one-sided formula for \code{svyqqmath} or a two-sided formula for \code{svyqqplot}
}
  \item{design}{
Survey design object to look up variables
}
  \item{designx}{
Survey design object to look up the RHS variable in \code{svyqqplot}, if
different from the LHS variable
}
\item{null}{Quantile function to compare the data quantiles to}
  \item{na.rm}{
Remove missing values
}
  \item{qrule}{
How to define quantiles for \code{svyqqplot} -- see
\code{\link{svyquantile}} for possible values
}
  \item{xlab,ylab}{
Passed to \code{plot}. For \code{svyqqplot}, if these are  \code{NULL}
they are replaced by the variable names
}

  \item{\dots}{
Graphical options to be passed to \code{plot}
}
}

\value{
None
}



\seealso{
\code{\link{quantile}}
\code{\link{qqnorm}}
\code{\link{qqplot}}
}
\examples{
data(api)

dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat,
fpc=~fpc)

svyqqmath(~api99, design=dstrat)
svyqqplot(api00~api99, design=dstrat)

dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
opar<-par(mfrow=c(1,2))

## sample distributions very different
qqplot(apiclus1$enroll, apistrat$enroll); abline(0,1)

## estimated population distributions much more similar
svyqqplot(enroll~enroll, design=dstrat,designx=dclus1,qrule=survey:::qrule_hf8); abline(0,1)
par(opar)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey }% use one of  RShowDoc("KEYWORDS")
