\name{crowd}
\alias{crowd}
\docType{data}
\title{Household crowding}
\description{
A tiny dataset from the VPLX manual.
}
\usage{data(crowd)}
\format{
  A data frame with 6 observations on the following 5 variables.
  \describe{
    \item{rooms}{Number of rooms in the house}
    \item{person}{Number of people in the household}
    \item{weight}{Sampling weight}
    \item{cluster}{Cluster number}
    \item{stratum}{Stratum number}
  }
}
\source{
  Manual for VPLX, Census Bureau.
}
\examples{
data(crowd)

## Example 1-1
i1.1<-as.svrepdesign(svydesign(id=~cluster, weight=~weight,data=crowd))
i1.1<-update(i1.1, room.ratio=rooms/person,
overcrowded=factor(person>rooms))
svymean(~rooms+person+room.ratio,i1.1)
svytotal(~rooms+person+room.ratio,i1.1)
svymean(~rooms+person+room.ratio,subset(i1.1,overcrowded==TRUE))
svytotal(~rooms+person+room.ratio,subset(i1.1,overcrowded==TRUE))

## Example 1-2
i1.2<-as.svrepdesign(svydesign(id=~cluster,weight=~weight,strata=~stratum, data=crowd))
svymean(~rooms+person,i1.2)
svytotal(~rooms+person,i1.2)

}
\keyword{datasets}
