\name{fpc}
\alias{fpc}
\non_function{}
\title{Small survey example}
\usage{data(fpc)}
\description{
The \code{fpc} data frame has 8 rows and 6 columns. It is artificial
data to illustrate survey sampling estimators.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{stratid}{Stratum ids}
    \item{psuid}{Sampling unit ids}
    \item{weight}{Sampling weights}
    \item{nh}{number sampled per stratum}
    \item{Nh}{population size per stratum}
    \item{x}{data}
  }
}

\source{
\url{http://www.stata-press.com/data/r7/fpc.dta}
}

\examples{
data(fpc)
fpc


withoutfpc<-svydesign(weights=~weight, ids=~psuid, strata=~stratid, variables=~x, data=fpc, nest=TRUE)

withoutfpc
svymean(~x, withoutfpc)

withfpc<-svydesign(weights=~weight, ids=~psuid, strata=~stratid,
fpc=~Nh, variables=~x, data=fpc, nest=TRUE)

withfpc
svymean(~x, withfpc)

## Other equivalent forms 
withfpc<-svydesign(prob=~I(1/weight), ids=~psuid, strata=~stratid,
fpc=~Nh, variables=~x, data=fpc, nest=TRUE)

svymean(~x, withfpc)

withfpc<-svydesign(weights=~weight, ids=~psuid, strata=~stratid,
fpc=~I(nh/Nh), variables=~x, data=fpc, nest=TRUE)

svymean(~x, withfpc)

withfpc<-svydesign(weights=~weight, ids=~interaction(stratid,psuid),
strata=~stratid, fpc=~I(nh/Nh), variables=~x, data=fpc)

svymean(~x, withfpc)

withfpc<-svydesign(ids=~psuid, strata=~stratid, fpc=~Nh,
 variables=~x,data=fpc,nest=TRUE)

svymean(~x, withfpc)

withfpc<-svydesign(ids=~psuid, strata=~stratid,
fpc=~I(nh/Nh), variables=~x, data=fpc, nest=TRUE)

svymean(~x, withfpc)



}
\keyword{datasets}
