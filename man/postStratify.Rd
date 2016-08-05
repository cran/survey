\name{postStratify}
\alias{postStratify}
\alias{postStratify.twophase}
\alias{postStratify.svyrep.design}
\alias{postStratify.survey.design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Post-stratify a survey }
\description{
 Post-stratification adjusts the sampling and replicate weights so that
 the joint distribution of a set of post-stratifying variables matches
 the known population joint distribution. Use \code{\link{rake}} when
 the full joint distribution is not available.
}
\usage{
postStratify(design, strata, population, partial = FALSE, ...)
\method{postStratify}{svyrep.design}(design, strata, population, partial = FALSE, compress=NULL,...)
\method{postStratify}{survey.design}(design, strata, population, partial = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{A survey design with replicate weights}
  \item{strata}{A formula or data frame of post-stratifying variables, which must not contain missing values. }
  \item{population}{A \code{\link{table}}, \code{\link{xtabs}} or \code{data.frame}
    with population frequencies }
  \item{partial}{if \code{TRUE}, ignore population strata not present in
    the sample}
  \item{compress}{Attempt to compress the replicate weight matrix? When
    \code{NULL} will attempt to compress if the original weight matrix
    was compressed}
  \item{...}{arguments for future expansion}
}
\details{
  The \code{population} totals can be specified as a table with the
  strata variables in the margins, or as a data frame where one column
  lists frequencies and the other columns list the unique combinations
  of strata variables (the format produced by \code{as.data.frame}
  acting on a \code{table} object).   A table must have named dimnames
  to indicate the variable names.
  
  Compressing the replicate weights will take time and may even
  increase memory use if there is actually little redundancy in the
  weight matrix (in particular if the post-stratification variables have
  many values and cut  across PSUs).
  
  If a \code{svydesign} object is to be converted to a replication
  design the post-stratification should be performed after conversion.
 
  
  The variance estimate for replication designs follows the same
  procedure as Valliant (1993) described for estimating totals. Rao et
  al (2002) describe this procedure for estimating functions (and also
  the GREG or g-calibration procedure, see \code{\link{calibrate}})
  
}

\note{
  If the sampling weights are already post-stratified there will be no
  change in point estimates after \code{postStratify} but the standard
  error estimates will decrease to correctly reflect the post-stratification.
}

\value{
 A new survey design object.
}
\references{
Valliant R (1993) Post-stratification and conditional variance
estimation. JASA 88: 89-96  
  
Rao JNK, Yung W, Hidiroglou MA (2002)   Estimating equations for the
analysis of survey data using poststratification information. Sankhya
64 Series A Part 2, 364-378.
}

\seealso{
  \code{\link{rake}}, \code{\link{calibrate}} for other things to do
  with auxiliary information
  
  \code{\link{compressWeights}} for information on compressing weights}
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1<-as.svrepdesign(dclus1)

svymean(~api00, rclus1)
svytotal(~enroll, rclus1)

# post-stratify on school type
pop.types <- data.frame(stype=c("E","H","M"), Freq=c(4421,755,1018))
#or: pop.types <- xtabs(~stype, data=apipop)
#or: pop.types <- table(stype=apipop$stype)

rclus1p<-postStratify(rclus1, ~stype, pop.types)
summary(rclus1p)
svymean(~api00, rclus1p)
svytotal(~enroll, rclus1p)


## and for svydesign objects
dclus1p<-postStratify(dclus1, ~stype, pop.types)
summary(dclus1p)
svymean(~api00, dclus1p)
svytotal(~enroll, dclus1p)
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}
