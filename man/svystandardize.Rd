\name{svystandardize}
\alias{svystandardize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Direct standardization within domains
}
\description{
In health surveys it is often of interest to standardize domains to have the same distribution of, eg, age as in a target population.  The operation is similar to post-stratification, except that the totals for the domains are fixed at the current estimates, not at known population values.  This function matches the estimates produced by the (US) National Center for Health Statistics.
}
\usage{
svystandardize(design, by, over, population, excluding.missing = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
 survey design object  
}
  \item{by}{
A one-sided formula specifying the variables whose distribution will be standardised
}
  \item{over}{
A one-sided formula specifying the domains within which the standardisation will occur
}
  \item{population}{
Desired population totals or proportions for the levels of combinations of variables in \code{by}
}
  \item{excluding.missing}{
  Optionally, a one-sided formula specifying variables whose missing values should be dropped before calculating the domain totals. 
}
}

\value{
A new survey design object of the same type as the input.
}
\references{
National Center for Health Statistics \url{http://www.cdc.gov/nchs/tutorials/NHANES/NHANESAnalyses/agestandardization/age_standardization_intro.htm}}
\note{
The standard error estimates do not exactly match the NCHS estimates
}



\seealso{
\code{\link{postStratify}}, \code{\link{svyby}}
}
\examples{
## matches http://www.cdc.gov/nchs/data/databriefs/db92_fig1.png
data(nhanes)
popage <- c( 55901 , 77670 , 72816 , 45364 )
design<-svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, data=nhanes, nest=TRUE)
stdes<-svystandardize(design, by=~agecat, over=~race+RIAGENDR, 
   population=popage, excluding.missing=~HI_CHOL)
svyby(~HI_CHOL, ~race+RIAGENDR, svymean, design=subset(stdes, agecat!="(0,19]"))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

