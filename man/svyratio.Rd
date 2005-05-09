\name{svyratio}
\alias{svyratio}
\alias{svrepratio}
\alias{print.svyratio}
\alias{print.svyratio_separate}
\alias{svyratio.survey.design}
\alias{svyratio.survey.design2}
\alias{svyratio.svyrep.design}
\alias{SE.svyratio}
\alias{predict.svyratio}
\alias{predict.svyratio_separate}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Ratio estimation}
\description{
Ratio estimation and estimates of totals based on ratios for complex
survey samples.
}
\usage{
\method{svyratio}{survey.design2}(numerator, denominator, design,separate=FALSE,...)
\method{svyratio}{svyrep.design}(numerator, denominator, design,...)
\method{predict}{svyratio}(object, total, se=TRUE,...)
\method{predict}{svyratio_separate}(object, total, se=TRUE,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{numerator}{formula, expression, or data frame giving numerator variable(s)}
  \item{denominator}{formula, expression, or data frame giving
    denominator variable(s) }
  \item{design}{from \code{svydesign} for \code{svyratio}, from
    \code{svrepdesign} for \code{svrepratio}}
  \item{object}{result of \code{svyratio}}
  \item{total}{vector of population totals for the denominator variables in
    \code{object}, or list of vectors of 
    population stratum totals if \code{separate=TRUE}}
  \item{se}{Return standard errors?}
  \item{separate}{Estimate ratio separately for strata}
  \item{...}{Other unused arguments for other methods}
}
\details{
  The separate ratio estimate of a total is the sum of ratio estimates
  in each stratum.  If the stratum totals supplied in the \code{total}
  argument and the strata in the design object both have names these
  names will be matched. If they do not have names it is important that
  the sample totals are supplied in the correct order, the same order
  as shown in the output of \code{summary(design)}.
  }
\value{
\code{svyratio} returns an object of class \code{svyratio}. The
\code{predict} method returns a matrix of population totals and
optionally a matrix of standard errors.
}
\references{Levy and Lemeshow. "Sampling of Populations" (3rd edition). Wiley}
\author{Thomas Lumley}

\seealso{\code{\link{svydesign}}
  
  \code{\link{svymean}} for estimating proportions and domain means
  
  \code{\link{calibrate}} for estimators related to the separate ratio estimator.
}

\examples{
data(scd)

## survey design objects
scddes<-svydesign(data=scd, prob=~1, id=~ambulance, strata=~ESA,
nest=TRUE, fpc=rep(5,6))
scdnofpc<-svydesign(data=scd, prob=~1, id=~ambulance, strata=~ESA,
nest=TRUE)

# convert to BRR replicate weights
scd2brr <- as.svrepdesign(scdnofpc, type="BRR")

# use BRR replicate weights from Levy and Lemeshow
repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
c(0,1,0,1,1,0))
scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)

# ratio estimates
svyratio(~alive, ~arrests, design=scddes)
svyratio(~alive, ~arrests, design=scdnofpc)
svyratio(~alive, ~arrests, design=scd2brr)
svyratio(~alive, ~arrests, design=scdrep)


data(api)
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)

## domain means are ratio estimates, but available directly
svyratio(~I(api.stu*(comp.imp=="Yes")), ~as.numeric(comp.imp=="Yes"), dstrat)
svymean(~api.stu, subset(dstrat, comp.imp=="Yes"))

## separate and combined ratio estimates of total
(sep<-svyratio(~api.stu,~enroll, dstrat,separate=TRUE))
(com<-svyratio(~api.stu, ~enroll, dstrat))

stratum.totals<-list(E=1877350, H=1013824, M=920298)

predict(sep, total=stratum.totals)
predict(com, total=sum(unlist(stratum.totals)))


}
\keyword{survey}% at least one, from doc/KEYWORDS

