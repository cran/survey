\name{svylogrank}
\alias{svylogrank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare survival distributions
}
\description{
Computes a weighted version of the logrank and stratified logrank tests
for comparing two or more survival distributions.  The generalization to
complex samples is based on the characterization of the logrank test as
the score test in a Cox model,  Under simple random sampling with
replacement, this function with \code{rho=0} and \code{gamma=0} 
is almost identical to the robust score test
in the survival package.
}
\usage{
svylogrank(formula, design, rho=0,gamma=0,method=c("small","large","score"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  Model formula with a single predictor and optionally a \code{\link{strata}} term. The predictor must be a factor if it has more than two levels.
}
  \item{design}{
  A survey design object
}
\item{rho,gamma}{Coefficients for the Harrington/Fleming G-rho-gamma
  tests. The default is the logrank test, \code{rho=1} gives a
  generalised Wilcoxon test}
\item{method}{\code{"small"} works faster when a matrix with dimension
  number of events by number of people fits easily in memory;
  \code{"large"} works faster for large data sets;  \code{"score"} works
by brute-force construction of an expanded data set, and is for debugging}
  \item{\dots}{
  for future expansion.
}
}

\value{
A vector containing the z-statistic for comparing each level of the variable to the lowest, the chisquared statistic for the logrank test, and the p-value.
}

\references{
Rader K (2014) 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{svykm}}, \code{\link{svycoxph}}.
}
\examples{
library("survival")
data(nwtco)
## stratified on case status
dcchs<-twophase(id=list(~seqno,~seqno), strata=list(NULL,~rel),
         subset=~I(in.subcohort | rel), data=nwtco, method="simple")
svylogrank(Surv(edrel,rel)~factor(stage),design=dcchs)

data(pbc, package="survival")
pbc$randomized <- with(pbc, !is.na(trt) & trt>0)
biasmodel<-glm(randomized~age*edema,data=pbc)
pbc$randprob<-fitted(biasmodel)
dpbc<-svydesign(id=~1, prob=~randprob, strata=~edema, data=subset(pbc,randomized))

svylogrank(Surv(time,status==2)~trt,design=dpbc)

svylogrank(Surv(time,status==2)~trt,design=dpbc,rho=1)

rpbc<-as.svrepdesign(dpbc)
svylogrank(Surv(time,status==2)~trt,design=rpbc)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{survival}% __ONLY ONE__ keyword per line
