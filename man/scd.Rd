\name{scd}
\alias{scd}
\non_function{}
\title{Survival in cardiac arrest}
\usage{data(scd)}
\description{
These data are from Section 12.2 of Levy and Lemeshow.  They describe
(a possibly apocryphal) study of survival in out-of-hospital cardiac
arrest. Two out of five ambulance stations were sampled from each of
three emergency service areas.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{ESA}{Emergency Service Area (strata)}
    \item{ambulance}{Ambulance station (PSU)}
    \item{arrests}{estimated number of cardiac arrests}
    \item{alive}{number reaching hospital alive}
  }
}

\source{
Levy and Lemeshow. "Sampling of Populations" (3rd edition). Wiley.
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
# or to Rao-Wu bootstrap
scd2boot <- as.svrepdesign(scdnofpc, type="subboot")

# use BRR replicate weights from Levy and Lemeshow
repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
c(0,1,0,1,1,0))
scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)

# ratio estimates
svyratio(~alive, ~arrests, design=scddes)
svyratio(~alive, ~arrests, design=scdnofpc)
svyratio(~alive, ~arrests, design=scd2brr)
svyratio(~alive, ~arrests, design=scd2boot)
svyratio(~alive, ~arrests, design=scdrep)

# or a logistic regression
summary(svyglm(cbind(alive,arrests-alive)~1, family=quasibinomial, design=scdnofpc))
summary(svyglm(cbind(alive,arrests-alive)~1, family=quasibinomial, design=scdrep))

# Because no sampling weights are given, can't compute design effects
# without replacement: use deff="replace"

svymean(~alive+arrests, scddes, deff=TRUE)
svymean(~alive+arrests, scddes, deff="replace")

}
\keyword{datasets}
