\name{svrepdesign}
\alias{svrepdesign}
\alias{svrepdesign.default}
\alias{svrepdesign.imputationList}
\alias{[.svyrep.design}
\alias{image.svyrep.design}
\alias{print.svyrep.design}
\alias{model.frame.svyrep.design}
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
svrepdesign(variables , repweights , weights, data,...)
\method{svrepdesign}{default}(variables = NULL, repweights = NULL, weights = NULL, data =
NULL, type = c("BRR", "Fay", "JK1","JKn","bootstrap","other"),
combined.weights=TRUE, rho = NULL, bootstrap.average=NULL,
scale=NULL, rscales=NULL,fpc=NULL, fpctype=c("fraction","correction"),...)
\method{svrepdesign}{imputationList}(variables=NULL, repweights,weights,data,...)
\method{image}{svyrep.design}(x, ..., col=grey(seq(.5,1,length=30)), type.=c("rep","total"))
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
  \item{bootstrap.average}{For \code{type="bootstrap"}, if the bootstrap
    weights have been averaged, gives the number of iterations averaged over}
  \item{scale, rscales}{Scaling constant for variance, see Details
below}
\item{fpc,fpctype}{Finite population correction information}
\item{x}{survey design with replicate weights}
\item{...}{Other arguments to \code{\link{image}}}
\item{col}{Colors}
\item{type.}{\code{"rep"} for only the replicate weights, \code{"total"} for the replicate and sampling weights combined.}
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

  Averaged bootstrap weights ("mean bootstrap") are used for some
  surveys from Statistics Canada. Yee et al (1999) describe their
  construction and use for one such survey.
  
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

  To generate your own replicate weights either use
  \code{\link{as.svrepdesign}} on a \code{survey.design} object, or see
  \code{\link{brrweights}}, \code{\link{bootweights}},
  \code{\link{jk1weights}} and \code{\link{jknweights}}
  
  The \code{model.frame} method extracts the observed data.

}
\value{
  Object of class \code{svyrep.design}, with methods for \code{print},
  \code{summary}, \code{weights}, \code{image}.
}
\references{Levy and Lemeshow. "Sampling of Populations". Wiley.

  Shao and Tu. "The Jackknife and Bootstrap." Springer.

  Yee et al (1999). Bootstrat Variance Estimation for the National
    Population Health Survey. Proceedings of the ASA Survey Research
    Methodology Section. \url{ http://www.amstat.org/Sections/Srms/Proceedings/papers/1999_136.pdf}
}

\note{To use replication-weight analyses on a survey specified by
  sampling design, use \code{as.svrepdesign} to convert it. }

\seealso{\code{\link{as.svrepdesign}}, \code{\link{svydesign}},
  \code{\link{brrweights}}, \code{bootweights} }

\examples{
data(scd)
# use BRR replicate weights from Levy and Lemeshow
repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
c(0,1,0,1,1,0))
scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights, combined.weights=FALSE)
svyratio(~alive, ~arrests, scdrep)
}
\keyword{survey}% at least one, from doc/KEYWORDS

