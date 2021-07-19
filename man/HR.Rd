\name{HR}
\Rdversion{1.1}
\alias{HR}
\alias{ppsmat}
\alias{ppscov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wrappers for specifying PPS designs}
\description{
The Horvitz-Thompson estimator and the Hartley-Rao approximation require information in addition to the sampling probabilities for sampled individuals.  These functions allow this information to be supplied.
}
\usage{
HR(psum=NULL, strata = NULL)
ppsmat(jointprob, tolerance = 1e-04)
ppscov(probcov, weighted=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psum}{ The sum of squared sampling probabilities for the population, divided by the sample size, as a single number or as a vector for stratified sampling
}
  \item{strata}{
Stratum labels, of the same length as \code{psum}, if \code{psum} is a vector
}
\item{jointprob}{Matrix of pairwise sampling probabilities for the sampled individuals}
\item{tolerance}{Tolerance for deciding that the covariance of sampling indicators is zero}
\item{probcov}{Covariance of the sampling indicators (often written 'Delta'), or weighted covariance if \code{weighted=TRUE}}
\item{weighted}{If \code{TRUE}, the \code{probcov} argument is the covariance divided by pairwise sampling probabilities}
}
\value{
An object of class \code{HR},\code{ppsmat}, \code{ppsdelta}, or \code{ppsdcheck} suitable for supplying as the \code{pps} argument to \code{\link{svydesign}}.
}

\seealso{
\link{election} for examples of PPS designs
}
\examples{
HR(0.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

