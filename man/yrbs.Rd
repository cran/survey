\name{yrbs}
\alias{yrbs}
\docType{data}
\title{
One variable from the Youth Risk Behaviors Survey, 2015.
}
\description{
Design information from the Youth Risk Behaviors Survey (YRBS), together
with the single variable `Never/Rarely wore bike helmet'.  Used as an
analysis example by CDC.
}
\usage{data("yrbs")}
\format{
  A data frame with 15624 observations on the following 4 variables.
  \describe{
    \item{\code{weight}}{sampling weights}
    \item{\code{stratum}}{sampling strata}
    \item{\code{psu}}{primary sampling units}
    \item{\code{qn8}}{1=Yes, 2=No}
  }
}

\source{
  \url{ftp://ftp.cdc.gov/pub/data/yrbs/2015smy/} for files
  
}
\references{
Centers for Disease Control and Prevention (2016) Software for Analysis
of YRBS Data. Previously at \url{https://www.cdc.gov/healthyyouth/data/yrbs/pdf/2015/2015_yrbs_analysis_software.pdf}
}
\examples{
data(yrbs)

yrbs_design <- svydesign(id=~psu, weight=~weight, strata=~stratum,
data=yrbs)
yrbs_design <- update(yrbs_design, qn8yes=2-qn8)

ci <- svyciprop(~qn8yes, yrbs_design, na.rm=TRUE, method="xlogit")
ci

## to print more digits: matches SUDAAN and SPSS exactly, per table 3 of reference
coef(ci)
SE(ci)
attr(ci,"ci")
}
\keyword{datasets}
