\name{postStratify}
\alias{postStratify}
\alias{postStratify.svyrep.design}
\alias{postStratify.survey.design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Post-stratify a survey }
\description{
 Post-stratification adjusts the sampling and replicate weights so that
 the joint distribution of a set of post-stratifying variables matches
 the known population joint distribution.
}
\usage{
postStratify(design, strata, population, partial = FALSE, ...)
\method{postStratify}{svyrep.design}(design, strata, population, partial = FALSE, compress=NULL,...)
\method{postStratify}{survey.design}(design, strata, population, partial = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{A survey design with replicate weights}
  \item{strata}{A formula or data frame of post-stratifying variables }
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
  
  The method for \code{survey.design} objects should be used only in
  fairly large samples as the resulting standard errors use a
  simple first-order approximation.
  
}
\value{
 A new survey design object.
}
\references{ }

\seealso{\code{\link{as.svrepdesign}}, \code{\link{svrepdesign}},
  \code{\link{rake}}, \code{\link{compressWeights}} }
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1<-as.svrepdesign(dclus1)

svrepmean(~api00, rclus1)
svreptotal(~enroll, rclus1)

# post-stratify on school type
pop.types <- data.frame(stype=c("E","H","M"), Freq=c(4421,755,1018))
#or: pop.types <- xtabs(~stype, data=apipop)
#or: pop.types <- table(stype=apipop$stype)

rclus1p<-postStratify(rclus1, ~stype, pop.types)
summary(rclus1p)
svrepmean(~api00, rclus1p)
svreptotal(~enroll, rclus1p)


## and for svydesign objects
dclus1p<-postStratify(dclus1, ~stype, pop.types)
summary(dclus1p)
svymean(~api00, dclus1p)
svytotal(~enroll, dclus1p)
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}
