\name{xdesign}
\alias{xdesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Crossed effects and other sparse correlations
}
\description{
  Defines a design object with multiple dimensions of correlation:
  observations that share any of the \code{id} variables are correlated,
  or you can supply an adjacency matrix or Matrix to specify which are
  correlated. Supports crossed designs (eg multiple raters of multiple
  objects) and non-nested observational correlation (eg observations
  sharing primary school or secondary school). Has methods for
  \code{svymean}, \code{svytotal}, \code{svyglm} (so far). 
}
\usage{
xdesign(id = NULL, strata = NULL, weights = NULL, data, fpc = NULL,
adjacency = NULL, overlap = c("unbiased", "positive"), allow.non.binary = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
list of formulas specifying cluster identifiers for each clustering dimension (or \code{NULL})
  }
  \item{strata}{
    Not implemented
  }
  \item{weights}{
    model formula specifying (sampling) weights
  }
  \item{data}{
    data frame containing all the variables
}
  \item{fpc}{
    Not implemented
  }
  \item{adjacency}{
    Adjacency matrix or Matrix indicating which pairs of observations
    are correlated
  }
  \item{overlap}{
   See details below
}
  \item{allow.non.binary}{

    If \code{FALSE} check that \code{adjacency} is a binary 0/1 or
    \code{TRUE}/\code{FALSE} matrix or Matrix.
}
}
\details{
  Subsetting for these objects actually drops observations; it is not
    equivalent to just setting weights to zero as for survey
    designs. So, for example, a subset of a balanced design will not be
    a balanced design.
  
  The \code{overlap} option controls double-counting of some variance
    terms. Suppose there are two clustering dimensions, \code{~a} and
    \code{~b}. If we compute variance matrices clustered on \code{a} and
    clustered on \code{b} and add them, observations that share both
    \code{a} and \code{b} will be counted twice, giving a positively
    biased estimator.  We can subtract off a variance matrix clustered
    on combinations of \code{a} and \code{b} to give an unbiased
    variance estimator.  However, the unbiased estimator is not
    guaranteed to be positive definite. In the references, Miglioretti
    and Heagerty use the \code{overlap="positive"} estimator and Cameron
    et al use the \code{overlap="unbiased"} estimator. 

  
}
\value{
An object of class \code{xdesign}
}
\references{
Miglioretti D, Heagerty PJ (2007) Marginal modeling of nonnested
multilevel data using standard software. Am J Epidemiol 165(4):453-63

Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2011). Robust Inference
With Multiway Clustering. Journal of Business & Economic Statistics,
29(2), 238-249.

\url{https://notstatschat.rbind.io/2021/09/18/crossed-clustering-and-parallel-invention/}
}

\seealso{\code{\link{salamander}}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{


## With one clustering dimension, is close to the with-replacement
##   survey estimator, but not identical unless clusters are equal size
data(api)
dclus1r<-svydesign(id=~dnum, weights=~pw, data=apiclus1)
xclus1<-xdesign(id=list(~dnum), weights=~pw, data=apiclus1)
xclus1

svymean(~enroll,dclus1r)
svymean(~enroll,xclus1)

data(salamander)
xsalamander<-xdesign(id=list(~Male, ~Female), data=salamander,
    overlap="unbiased")
xsalamander
degf(xsalamander)

}