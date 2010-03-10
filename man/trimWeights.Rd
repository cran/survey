\name{trimWeights}
\Rdversion{1.1}
\alias{trimWeights}
\alias{trimWeights.svyrep.design}
\alias{trimWeights.survey.design2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trim sampling weights
}
\description{
Trims very high or very low sampling weights to reduce the influence of outlying observations. In a replicate-weight design object, the replicate weights are also trimmed. The total amount trimmed is divided among the observations that were not trimmed, so that the total weight remains the same.
}
\usage{
trimWeights(design, upper = Inf, lower = -Inf, ...)
\method{trimWeights}{survey.design2}(design, upper = Inf, lower = -Inf, strict=FALSE,...)
\method{trimWeights}{svyrep.design}(design, upper = Inf, lower = -Inf,compress=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
  A survey design object
}
  \item{upper}{
  Upper bound for weights
}
  \item{lower}{
  Lower bound for weights
  }
  \item{strict}{
    The reapportionment of the `trimmings' from the weights can push
    other weights over the limits. If \code{trim=TRUE} the function
    calls itself recursively to prevent this.
  }
  \item{compress}{
    Compress the replicate weights after trimming.
    }
  \item{\dots}{
  Other arguments for future expansion
}
}

\value{
A new survey design object with trimmed weights.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{calibrate}} has a \code{trim} option for trimming the
calibration adjustments.
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

