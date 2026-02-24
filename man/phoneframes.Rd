\docType{data}
\name{phoneframes}
\alias{phoneframes}
\alias{DatA}
\alias{DatB}
\alias{PiklA}
\alias{PiklB}
\title{Database of household expenses for two sampling frames}
\format{\describe{
   \item{Domain}{A factor indicating the domain each household belongs
     to. In sample A, possible values are "a" if household belongs to
     domain a or "ab" if household belongs to overlap domain; in sample
     B, the values are "b" or "ba"}
   \item{Feed}{Feeding expenses (in euros) at the househould}
   \item{Clo}{Clothing expenses (in euros) at the household}
   \item{Lei}{Leisure expenses (in euros) at the household}
   \item{Inc}{Household income (in euros). Values for this variable are only available for households included in frame A. For households included in domain b, value of this variable is missing}
   \item{Tax}{Household municipal taxes (in euros) paid. Values for this
     variable are only available for households included in frame A. For
     households included in domain b, value of this variable is missing}
   \item{M2}{Square meters of the house. Values for this variable are only available for households included in frame B. For households included in domain a, value of this variable is missing}
   \item{Size}{Household size. Values for this variable are only available for households included in frame B. For households included in domain a, value of this variable is missing}
   \item{ProbA}{First order inclusion probability in frame A. This probability is 0 for households included in domain b.}
   \item{ProbB}{First order inclusion probability in frame B. This probability is 0 for households included in domain a.}
   \item{Stratum}{A numeric value indicating the stratum each household belongs to.}
}}
\usage{
data(phoneframes)
}
\description{
This dataset contains some variables regarding household expenses for a
sample of 105 households selected from a list of landline phones (frame
A) and a sample of 135 from a list of mobile phones (frame B) in a
particular city in a specific month. These data are taken from the
\code{Frames2} package under the GPL-2 or GPL-3 licence. 
}
\details{
The frame A sample, of size \eqn{n_A = 105}, has been drawn from a population of \eqn{N_A = 1735} households with landline phone according to a stratified random sampling. Population units were divided in 6 different strata. 
 Population sizes of these strata are \eqn{N_A^h = (727, 375, 113, 186,
   115, 219)}. \eqn{N_{ab} = 601} of the households composing the
 population have, also, mobile phone. On the other hand, frame totals
 for auxiliary variables in this frame are \eqn{X_{Income}^A = 4300260}
 and \eqn{X_{Taxes}^A = 215577}.

 The frame B sample, of size \eqn{n_B = 135}, has been drawn from a population of \eqn{N_B = 1191} households with mobile phone according to a simple random sampling without replacement design. 
 \eqn{N_{ab} = 601} of these households have, also, landline phone. On the other hand, frame totals for auxiliary variables in this frame are \eqn{X_{Metres2}^B = 176553} and \eqn{X_{Size}^B = 3529}

 \code{PiklA} and \code{PiklB} are matrices of pairwise sampling
 probabilities for the two frames. 
}
\examples{
data(phoneframes)
A_in_frames<-cbind(1, DatA$Domain=="ab")
B_in_frames<-cbind(DatB$Domain=="ba",1)

Bdes_pps<-svydesign(id=~1, fpc=~ProbB, data=DatB,pps=ppsmat(PiklB))
Ades_pps <-svydesign(id=~1, fpc=~ProbA,data=DatA,pps=ppsmat(PiklA))

## optimal constant (Hartley) weighting
mf_pps<-multiframe(list(Ades_pps,Bdes_pps),list(A_in_frames,B_in_frames),theta=0.74) 
svytotal(~Lei,mf_pps)

Awts<-cbind(1/DatA$ProbA, ifelse(DatA$ProbB==0,0,1/DatA$ProbB))
Bwts<-cbind(ifelse(DatB$ProbA==0,0,1/DatB$ProbA),1/DatB$ProbB )
## dividing by the expected number of selections (BKA or HH estimator)
mf_pps2<-multiframe(list(Ades_pps,Bdes_pps),list(Awts,Bwts),estimator="expected") 
svymean(~Lei,mf_pps2)

## Metcalf and Scott approximation
DatB$Stratum<-10
DatB$Frame<-2
DatA$Frame<-1
Dat_both<-rbind(DatA,DatB)
frame_weights<-c(0.742,1-0.742)
Dat_both$fweights<-with(Dat_both, ifelse(Frame==1, 
   ifelse(Domain=="ab", frame_weights[1]*1/ProbA,1/ProbA),
   ifelse(Domain=="ba", frame_weights[2]*1/ProbB, 1/ProbB)))

MSdesign<-svydesign(id=~1, strata=~Stratum, weights=~fweights,data=Dat_both)
svymean(~Lei,MSdesign)

}
\seealso{
  \code{\link{multiframe}}, \code{\link{reweight}}

Original package: \url{https://CRAN.R-project.org/package=Frames2}
}
\keyword{datasets}

