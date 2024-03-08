\name{svyrecvar}
\alias{svyrecvar}
\alias{multistage}
\alias{multistage_rcpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variance estimation for multistage surveys}
\description{
  Compute the variance of a total under multistage sampling, using a
  recursive descent algorithm.  
}
\usage{
svyrecvar(x, clusters, stratas,fpcs, postStrata = NULL,
lonely.psu = getOption("survey.lonely.psu"),
one.stage=getOption("survey.ultimate.cluster"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Matrix of data or estimating functions}
  \item{clusters}{Data frame or matrix with cluster ids for each stage}
  \item{stratas}{Strata for each stage }
  \item{fpcs}{Information on population and sample size for each stage,
    created by \code{\link{as.fpc}}}
  \item{postStrata}{post-stratification information as created by
    \code{\link{postStratify}} or \code{\link{calibrate}} }
  \item{lonely.psu}{How to handle strata with a single PSU}
  \item{one.stage}{If \code{TRUE}, compute a one-stage
    (ultimate-cluster) estimator}
}
\details{
  The main use of this function is to compute the variance of the sum
  of a set of estimating functions under multistage sampling.  The
  sampling is assumed to be simple or stratified random sampling within
  clusters at each stage except perhaps the last stage.  The variance of
  a statistic is computed from the variance of estimating functions as
  described by Binder (1983).
  
  Use \code{one.stage=FALSE} for compatibility with other software that
  does not perform multi-stage calculations, and set
  \code{options(survey.ultimate.cluster=TRUE)} to make this the default.
  
  The idea of a recursive algorithm is due to Bellhouse (1985).
  Texts such as Cochran (1977) and Sarndal et al (1991) describe the
  decomposition of the variance into a single-stage between-cluster
  estimator and a within-cluster estimator, and this is applied recursively.
  
  If \code{one.stage} is a positive integer it specifies the number of
  stages of sampling to use in the recursive estimator.

  If \code{pps="brewer"}, standard errors are estimated using Brewer's
  approximation for PPS without replacement, option 2 of those described
  by Berger (2004). The \code{fpc} argument must then be specified in
  terms of sampling fractions, not population sizes (or omitted, but
  then the \code{pps} argument would have no effect and the
  with-replacement standard errors would be correct).
}
\value{
   A covariance matrix
}
\references{
  Bellhouse DR (1985) Computing Methods for Variance Estimation in Complex Surveys.
  Journal of Official Statistics. Vol.1, No.3, 1985

  Berger, Y.G. (2004), A Simple Variance Estimator for Unequal
  Probability Sampling Without Replacement. Journal of Applied
  Statistics, 31, 305-315.
  
  Binder, David A. (1983).  On the variances of asymptotically normal
  estimators from complex surveys.  International Statistical Review,
  51, 279-292.

  Brewer KRW (2002) Combined Survey Sampling Inference (Weighing Basu's
  Elephants)  [Chapter 9]

  Cochran, W. (1977)  Sampling Techniques. 3rd edition. Wiley.
  
  Sarndal C-E, Swensson B, Wretman J (1991) Model Assisted Survey
  Sampling. Springer.
  
}
\note{
  A simple set of finite population corrections will only be exactly
  correct when each successive stage uses simple or stratified random
  sampling without replacement.  A correction under general unequal
  probability sampling (eg PPS) would require joint inclusion probabilities (or,
  at least, sampling probabilities for units not included in the sample),
  information not generally available.

  The quality of Brewer's approximation is excellent in Berger's
  simulations, but the accuracy may vary depending on the sampling
  algorithm used.
}


\seealso{
  \code{\link{svrVar}} for replicate weight designs
    
  \code{\link{svyCprod}} for a description of how variances are
    estimated at each stage
  }
  \examples{
data(mu284)
dmu284<-svydesign(id=~id1+id2,fpc=~n1+n2, data=mu284)
svytotal(~y1, dmu284)


data(api)
# two-stage cluster sample
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)
summary(dclus2)
svymean(~api00, dclus2)
svytotal(~enroll, dclus2,na.rm=TRUE)

# bootstrap for multistage sample
mrbclus2<-as.svrepdesign(dclus2, type="mrb", replicates=100)
svytotal(~enroll, mrbclus2, na.rm=TRUE)

# two-stage `with replacement'
dclus2wr<-svydesign(id=~dnum+snum, weights=~pw, data=apiclus2)
summary(dclus2wr)
svymean(~api00, dclus2wr)
svytotal(~enroll, dclus2wr,na.rm=TRUE)


}
\keyword{survey}% at least one, from doc/KEYWORDS

