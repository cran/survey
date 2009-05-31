\name{election}
\alias{election}
\alias{election_pps}
\alias{election_jointprob}
\alias{election_jointHR}
\alias{election_insample}
\docType{data}
\title{US 2004 presidential election data at state or county level}
\description{
A sample of voting data from US states or counties (depending on data
availability), sampled with probability proportional to number of votes. The sample was drawn using Tille's splitting method, implemented in the "sampling" package.
}
\usage{data(election)}
\format{
  \code{election} is a data frame with 4600 observations on the following 8 variables.
  \describe{
    \item{\code{County}}{A factor specifying the state or country}
    \item{\code{TotPrecincts}}{Number of precincts in the state or county}
    \item{\code{PrecinctsReporting}}{Number of precincts supplying data}
    \item{\code{Bush}}{Votes for George W. Bush}
    \item{\code{Kerry}}{Votes for John Kerry}
    \item{\code{Nader}}{Votes for Ralph Nader}
    \item{\code{votes}}{Total votes for those three candidates}
    \item{\code{p}}{Sampling probability, proportional to \code{votes}}
  }
  
  \code{election_pps} is a sample of 40 counties or states taken with
  probability proportional to the number of votes. It includes the
  additional column \code{wt} with the sampling weights. 

  \code{election_insample} indicates which rows of \code{election} were sampled.

  \code{election_jointprob} are the pairwise sampling probabilities and
  \code{election_jointHR} are approximate pairwise sampling probabilities using
  the Hartley-Rao approximation.
}
\source{
.
}

\examples{
data(election)
## high positive correlation between totals
plot(Bush~Kerry,data=election,log="xy")
## high negative correlation between proportions
plot(I(Bush/votes)~I(Kerry/votes), data=election)

## Variances without replacement
## Horvitz-Thompson type
dpps_br<- svydesign(id=~1,  fpc=~p, data=election_pps, pps="brewer")
dpps_ov<- svydesign(id=~1,  fpc=~p, data=election_pps, pps="overton")
dpps_hr<- svydesign(id=~1,  fpc=~p, data=election_pps, pps=HR(sum(election$p^2)/40))
dpps_hr1<- svydesign(id=~1, fpc=~p, data=election_pps, pps=HR())
dpps_ht<- svydesign(id=~1,  fpc=~p, data=election_pps, pps=ppsmat(election_jointprob))
## Yates-Grundy type
dpps_yg<- svydesign(id=~1,  fpc=~p, data=election_pps, pps=ppsmat(election_jointprob),variance="YG")
dpps_hryg<- svydesign(id=~1,  fpc=~p, data=election_pps, pps=HR(sum(election$p^2)/40),variance="YG")

## The with-replacement approximation
dppswr <-svydesign(id=~1, probs=~p, data=election_pps)

svytotal(~Bush+Kerry+Nader, dpps_ht)
svytotal(~Bush+Kerry+Nader, dpps_yg)
svytotal(~Bush+Kerry+Nader, dpps_hr)
svytotal(~Bush+Kerry+Nader, dpps_hryg)
svytotal(~Bush+Kerry+Nader, dpps_hr1)
svytotal(~Bush+Kerry+Nader, dpps_br)
svytotal(~Bush+Kerry+Nader, dpps_ov)
svytotal(~Bush+Kerry+Nader, dppswr)
}
\keyword{datasets}
\keyword{survey}
