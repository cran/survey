\name{multiphase}
\alias{multiphase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiphase sampling designs
}
\description{
These objects represent designs with arbitrarily many nested phases of
sampling, allowing estimation and calibration/raking at
each phase
}
\usage{
multiphase(ids, subset, strata, probs, data, fpc = NULL,
check.variable.phase=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ids}{
List of as many model formulas as phases describing ids for each phase. Each formula may
indicate multistage sampling
  }
  \item{subset}{
list of model formulas for each phase except the first, specifying a
logical vector of which observations from the previous phase are included
  }
  \item{strata}{
List of as many model formulas as phases describing strata for each phase. Each formula may
indicate multistage sampling, or \code{NULL} for no strata
  }
  \item{probs}{
List of as many model formulas or \code{pps_spec} objects as phases
describing sampling probabilities for each phase. Each formula may
indicate multistage sampling. Typically will either be \code{NULL} except for
phase 1 if \code{strata} are specified, or a matrix of class
\code{pps_spec} specifying pairwise probabilities or covariances. Use
\code{~1} at phase 1 to specify iid sampling from a generating model.
}
\item{data}{
data frame of data}
  \item{fpc}{
Finite population correction for the first phase, if needed
  }
\item{check.variable.phase}{Work out which phase each variable is
  observed in by looking at missing value patterns. You may want
  \code{FALSE} for simulations where the values aren't actually missing}
}
\details{
  Variance calculation uses a decomposition with sampling contributions
  at each stage, which are returned as the \code{phases} attribution of
  a variance-covariance matrix. The computations broadly follow the
  description for two-phase sampling in chapter 9 of Sarndal et al
  (1991); there is more detail in the vignette
}
\value{
Object of class \code{multiphase}
}
\references{
Sarndal, Swensson, and Wretman (1991) "Model Assisted Survey Sampling"
(Chapter 9)
}

\note{
There are currently methods for \code{svytotal},
\code{svymean}, \code{svyglm}, \code{svyvar}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{twophase} for older implementations of two-phase sampling

  \code{vignette("multiphase")} for computational details
  
}
\examples{
data(nwtco)
dcchs<-twophase(id = list(~seqno, ~seqno), strata = list(NULL, ~rel), 
    subset = ~I(in.subcohort | rel), data = nwtco)
mcchs<-multiphase(id = list(~seqno, ~seqno), strata = list(NULL, ~rel), 
    subset = list(~I(in.subcohort | rel)), probs = list(~1, NULL), 
    data = nwtco)
dcchs
mcchs
svymean(~edrel, dcchs)
svymean(~edrel, mcchs)

summary(svyglm(edrel~rel+histol+stage, design=dcchs))
summary(svyglm(edrel~rel+histol+stage, design=mcchs))

m<-calibrate(mcchs,~factor(stage)+rel, phase=2, calfun="raking")
vcov(svytotal(~factor(stage), m))
}