\name{nhanes}
\alias{nhanes}
\docType{data}
\title{
Cholesterol data from a US survey
}
\description{
Data extracted from NHANES 2009-2010 on high cholesterol.  
}
\usage{data(nhanes)}
\format{
  A data frame with 8591 observations on the following 7 variables.
  \describe{
    \item{\code{SDMVPSU}}{Primary sampling units}
    \item{\code{SDMVSTRA}}{Sampling strata}
    \item{\code{WTMEC2YR}}{Sampling weights}
    \item{\code{HI_CHOL}}{Numeric vector: 1 for total cholesterol over
      240mg/dl, 0 under 240mg/dl}
    \item{\code{race}}{1=Hispanic, 2=non-Hispanic white, 3=non-Hispanic
      black, 4=other}
    \item{\code{agecat}}{Age group \code{(0,19]} \code{(19,39]} \code{(39,59]} \code{(59,Inf]}}
    \item{\code{RIAGENDR}}{Gender: 1=male, 2=female}
  }
}
\source{
  \url{https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=laboratory&CycleBeginYear=2009}
}
\examples{
data(nhanes)
design <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE,data=nhanes)
design
}
\keyword{datasets}
