\name{mu284}
\alias{mu284}
\docType{data}
\title{Two-stage sample from MU284}
\description{
The MU284 population comes from Sarndal et al, and the complete data are
available from Statlib. These data are a two-stage sample from the
population, analyzed on page 143 of the book.
}
\usage{data(mu284)}
\format{
  A data frame with 15 observations on the following 5 variables.
  \describe{
    \item{\code{id1}}{identifier for PSU}
    \item{\code{n1}}{number of PSUs in population}
    \item{\code{id2}}{identifier for second-stage unit}
    \item{\code{y1}}{variable to be analysed}
    \item{\code{n2}}{number of second-stage units in this PSU}
  }
}
\source{
  Carl Erik Sarndal, Bengt Swensson, Jan Wretman. (1991) "Model Assisted
  Survey Sampling" Springer.
}
\references{
Full MU284 population at \url{http://lib.stat.cmu.edu/datasets/mu284}
}
\examples{
data(mu284)
(dmu284<-svydesign(id=~id1+id2,fpc=~n1+n2, data=mu284))
(ytotal<-svytotal(~y1, dmu284))
vcov(ytotal)
}
\keyword{datasets}
