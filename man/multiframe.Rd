\name{multiframe}
\alias{multiframe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Dual-frame and multi-frame surveys
}
\description{
Given a list of samples from K different sampling frames and information
about which observations are in which frame, constructs an object
representing the whole multi-frame sample. If an unit is in the
overlap of multiple frames in the population it is effectively split
into multiple separate units and so the weight is split if it is
sampled. To optimise the split of frame weights, see \code{\link{reweight}}
}
\usage{
multiframe(designs, overlaps, estimator = c("constant", "expected"), theta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{designs}{
List of survey design objects
  }
  \item{overlaps}{
list of matrices. Each matrix has K columns indicating whether the
observation is in frames 1-K. For the 'constant'-type estimator, this is
binary. For the \code{expected} estimator the entry in row i
and column k is the weight or probability that observation i would have had if sampled
from frame k. (weights if all >=1, probabilities if all <=1)
  }
  \item{estimator}{
    \code{"constant"} specifies Hartley's estimators in which the
    partition of weights is the same for each observation.
    \code{"expected"} weights each observation by the reciprocal of the
    expected number of times it is sampled; it is the estimator proposed
    by Bankier and by Kalton and Anderson. 
  }
  \item{theta}{
  Scale factors adding to 1 for splitting the overlap between frames
  }
}
\details{
It is not necessary that the frame samples contain exactly the same
variables or that they are in the same order, although only variables present in all the samples can be
used.  It is important that factor variables existing across more than
one frame sample have the same factor levels in all the samples. 
  
All these estimators assume sampling is independent between frames, and that
any observation sampled more than once is present in the dataset each time
it is sampled.
  }
\value{
Object of class \code{multiframe}
}
\references{
Bankier, M. D. (1986) Estimators Based on Several Stratified Samples With Applications to Multiple Frame Surveys. Journal of the American Statistical Association, Vol. 81, 1074 - 1079.

Hartley, H. O. (1962) Multiple Frames Surveys. Proceedings of the American Statistical Association, Social Statistics Sections, 203 - 206.
Hartley, H. O. (1974) Multiple frame methodology and selected applications. Sankhya C, Vol. 36, 99 - 118.

Kalton, G. and Anderson, D. W. (1986) Sampling Rare Populations. Journal of the Royal Statistical Society, Ser. A, Vol. 149, 65 - 82.

}

\seealso{
  \code{\link{svydesign}}, \code{\link{reweight}}

  For a simple introduction:
  Metcalf P and Scott AJ (2009) Using multiple frames in health
  surveys. Stat Med 28:1512-1523

  For general reference:
  Lohr SL, Rao JNK. Inference from dual frame surveys. Journal of the American Statistical Association 2000; 94:271-280.

  Lohr SL, Rao JNK. Estimation in multiple frame surveys. Journal of the American Statistical Association 2006; 101:1019-1030.
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
svymean(~Lei, mf_pps)

svyby(~Lei, ~Size, svymean, design=mf_pps)
svytable(~Size+I(Lei>20), mf_pps,round=TRUE)

}
