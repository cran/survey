\name{svrepdesign}
\alias{svrepdesign}
\alias{[.svyrep.design}
\alias{image.svyrep.design}
\alias{weights.svyrep.design}
\alias{print.svyrep.design}
\alias{summary.svyrep.design}
\alias{print.summary.svyrep.design}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Specify survey design with replicate weights}
\description{
Some recent large-scale surveys specify replication weights rather than
the sampling design (partly for privacy reasons).  This function specifies the
data structure for such a survey. 
}
\usage{
svrepdesign(variables = NULL, repweights = NULL, weights = NULL, data =
NULL, type = c("BRR", "Fay", "JK1","JKn","other"),
combined.weights=FALSE, rho = NULL, 
scale=NULL, rscales=NULL,fpc=NULL, fpctype=c("fraction","correction"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{variables}{formula or data frame specifying variables to include in the design (default is all) }
  \item{repweights}{formula or data frame specifying replication weights }
  \item{weights}{sampling weights }
  \item{data}{data frame to look up variables in formulas}
  \item{type}{Type of replication weights}
  \item{combined.weights}{\code{TRUE} if the \code{repweights} already
    include the sampling weights}
  \item{rho}{Shrinkage factor for weights in Fay's method}
  \item{scale, rscales}{Scaling constant for variance, see Details
below}
\item{fpc,fpctype}{Finite population correction information}
}
\details{
  In the BRR method, the dataset is split into halves, and the
  difference between halves is used to estimate the variance. In Fay's
  method, rather than removing observations from half the sample they
  are given weight \code{rho} in one half-sample and \code{2-rho} in the
  other.  The ideal BRR analysis is restricted to a design where each
  stratum has two PSUs, however, it has been used in a much wider class
  of surveys.
  
  The JK1 and JKn types are both jackknife estimators deleting one
  cluster at a time. JKn is designed for stratified and JK1 for
  unstratified designs. 
  
  The variance is computed as the sum of squared deviations of the
  replicates from their mean.  This may be rescaled: \code{scale} is an
  overall multiplier and \code{rscale} is a vector of
  replicate-specific multipliers for the squared deviations.  If the
  replication weights incorporate the sampling weights
  (\code{combined.weights=TRUE}) or for \code{type="other"} these must
  be specified, otherwise they can be guessed from the weights.


A finite population correction may be specified for \code{type="other"},
\code{type="JK1"} and \code{type="JKn"}.  \code{fpc} must be a vector
with one entry for each replicate. To specify sampling fractions use
\code{fpctype="fraction"} and to specify the correction directly use
\code{fpctype="correction"}


}
\value{
  Object of class \code{svyrep.design}, with methods for \code{print},
  \code{summary}, \code{weights}, \code{image}.
}
\references{Levy and Lemeshow. "Sampling of Populations". Wiley. }
\note{To use replication-weight analyses on a survey specified by
  sampling design, use \code{as.svrepdesign} to convert it. }

\seealso{\code{\link{as.svrepdesign}}, \code{\link{svydesign}}, \code{\link{brrweights}} }

\examples{

}
\keyword{survey}% at least one, from doc/KEYWORDS

