\name{withReplicates}
\alias{withReplicates}
\alias{withReplicates.svyrep.design}
\alias{withReplicates.svrepvar}
\alias{withReplicates.svrepstat}
\alias{vcov.svyrep.design}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Compute variances by replicate weighting}
\description{
Given a function or expression computing a statistic based on sampling
weights, \code{withReplicates} evaluates the statistic and produces a
replicate-based estimate of variance. \code{vcov.svrep.design} produces
the variance estimate from a set of replicates and the design object.
}
\usage{
withReplicates(design, theta,..., return.replicates=FALSE)
\method{withReplicates}{svyrep.design}(design, theta, rho = NULL, ..., 
     scale.weights=FALSE, return.replicates=FALSE)
\method{withReplicates}{svrepvar}(design, theta,  ...,  return.replicates=FALSE)
\method{withReplicates}{svrepstat}(design, theta,  ...,  return.replicates=FALSE)
\method{vcov}{svyrep.design}(object, replicates, centre,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{design}{A survey design with replicate weights (eg from \code{\link{svrepdesign}}) or a suitable object with replicate parameter estimates}
  \item{theta}{A function or expression: see Details below}
  \item{rho}{If \code{design} uses BRR weights, \code{rho} optionally
    specifies the parameter for Fay's variance estimator.}
  \item{\dots}{Other arguments to \code{theta}}
  \item{scale.weights}{Divide the probability weights by their sum (can
    help with overflow problems)}
  \item{return.replicates}{Return the replicate estimates as well as
    the variance?}
    \item{object}{The replicate-weights design object used to create the replicates}
    \item{replicates}{A set of replicates}
    \item{centre}{The centering value for variance calculation. If
      \code{object$mse} is \code{TRUE} 
      this is the result of estimation using the sampling weights, and
      must be supplied.  If  \code{object$mse} is \code{FALSE} the
      mean of the replicates is used and this argument is silently ignored.
}}
\details{
The method for \code{svyrep.design} objects evaluates a function or
expression using the sampling weights and then each set of replicate
weights.  The method for \code{svrepvar} objects evaluates the function
or expression on an estimated population covariance matrix and its
replicates, to simplify multivariate statistics such as structural
equation models.

For the \code{svyrep.design} method, if \code{theta} is a function its first argument will be a vector of
  weights and the second argument will be a data frame containing the
  variables from the design object.   If it is an expression, the sampling weights will be available as the
  variable \code{.weights}.  Variables in the design object will also
  be in scope.  It is possible to use global variables in the
  expression, but unwise, as they may be masked by local variables
  inside \code{withReplicates}.

For the \code{svrepvar} method a function will get the covariance
matrix as its first argument, and an expression will be evaluated with
\code{.replicate} set to the variance matrix.

For the \code{svrepstat} method a function will get the point estimate, and an expression will be evaluated with
\code{.replicate} set to each replicate.  The method can only be used
when the \code{svrepstat} object includes replicates.

}
\value{
  If \code{return.replicates=FALSE},  the weighted statistic, with the
  variance matrix as the \code{"var"} attribute. If
  \code{return.replicates=TRUE}, a list with elements \code{theta} for
  the usual return value and \code{replicates} for the replicates.
}
\seealso{ \code{\link{svrepdesign}}, \code{\link{as.svrepdesign}}, \code{\link{svrVar}}}

\examples{
data(scd)
repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
c(0,1,0,1,1,0))
scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)

a<-svyratio(~alive, ~arrests, design=scdrep)
print(a$ratio)
print(a$var)
withReplicates(scdrep, quote(sum(.weights*alive)/sum(.weights*arrests)))
withReplicates(scdrep, function(w,data)
sum(w*data$alive)/sum(w*data$arrests))

data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1<-as.svrepdesign(dclus1)
varmat<-svyvar(~api00+api99+ell+meals+hsg+mobility,rclus1,return.replicates=TRUE)
withReplicates(varmat, quote( factanal(covmat=.replicate, factors=2)$unique) )


data(nhanes)
nhanesdesign <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE,data=nhanes)
logistic <- svyglm(HI_CHOL~race+agecat+RIAGENDR, design=as.svrepdesign(nhanesdesign),
family=quasibinomial, return.replicates=TRUE)
fitted<-predict(logistic, return.replicates=TRUE, type="response")
sensitivity<-function(pred,actual) mean(pred>0.1 & actual)/mean(actual)
withReplicates(fitted, sensitivity, actual=logistic$y)

\dontrun{
library(quantreg)
data(api)
## one-stage cluster sample
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
## convert to bootstrap
bclus1<-as.svrepdesign(dclus1,type="bootstrap", replicates=100)

## median regression
withReplicates(bclus1, quote(coef(rq(api00~api99, tau=0.5, weights=.weights))))
}


## pearson correlation
dstrat <- svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
bstrat<- as.svrepdesign(dstrat,type="subbootstrap")

v <- svyvar(~api00+api99, bstrat, return.replicates=TRUE)
vcor<-cov2cor(as.matrix(v))[2,1]
vreps<-v$replicates
correps<-apply(vreps,1, function(v) v[2]/sqrt(v[1]*v[4]))

vcov(bstrat,correps, centre=vcor)


}
\keyword{survey}% at least one, from doc/KEYWORDS


