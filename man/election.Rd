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
availability), sampled with probability proportional to number of votes
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
  additional column \code{wt} with the sampling weights

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

dpps<- svydesign(id=~1, weights=~wt, fpc=~p, data=election_pps, pps="brewer")
dppswr <-svydesign(id=~1, weights=~wt, data=election_pps)
svytotal(~Bush+Kerry+Nader, dpps)
svytotal(~Bush+Kerry+Nader, dppswr)
}
\keyword{datasets}
\keyword{survey}
