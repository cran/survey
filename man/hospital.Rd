\name{hospital}
\alias{hospital}
\non_function{}
\title{Sample of obstetric hospitals }
\usage{data(hospital)}
\description{
The \code{hospital} data frame has 15 rows and 5 columns. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{hospno}{a numeric vector}
    \item{oblevel}{a numeric vector}
    \item{weighta}{a numeric vector}
    \item{tothosp}{a numeric vector}
    \item{births}{a numeric vector}
  }
}
\source{
\url{http://www.ats.ucla.edu/stat/books/sop/hospsamp.dta}
}
\references{
Levy and Lemeshow. "Sampling of Populations" (3rd edition). Wiley.
}
\examples{
data(hospital)
hospdes<-svydesign(strata=~oblevel, id=~hospno, weights=~weighta,
fpc=~tothosp, data=hospital)
hosprep<-as.svrepdesign(hospdes)

svytotal(~births, design=hospdes)
svreptotal(~births, design=hosprep)
}
\keyword{datasets}
