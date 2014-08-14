\name{as.svrepdesign}
\alias{as.svrepdesign}
%- Also NEED an `\alias' for EACH other topic documented here. 
\title{Convert a survey design to use replicate weights}
\description{
  Creates a replicate-weights survey design object from a traditional
  strata/cluster survey design object. \code{JK1} and \code{JKn} are
  jackknife methods, \code{BRR} is Balanced Repeated Replicates and
  \code{Fay} is Fay's modification of this, \code{bootstrap} is Canty
  and Davison's bootstrap, \code{subbootstrap} is Rao and Wu's
  \eqn{(n-1)} bootstrap, and \code{mrbbootstrap} is Preston's multistage rescaled bootstrap.
}
\usage{
as.svrepdesign(design, type=c("auto", "JK1", "JKn", "BRR", "bootstrap",
   "subbootstrap","mrbbootstrap","Fay"),
   fay.rho = 0, fpc=NULL,fpctype=NULL,..., compress=TRUE, 
   mse=getOption("survey.replicates.mse"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{design}{Object of class \code{survey.design} }
  \item{type}{Type of replicate weights. \code{"auto"} uses JKn for
    stratified, JK1 for unstratified designs}
  \item{fay.rho}{Tuning parameter for Fay's variance method }
  \item{fpc,fpctype,\dots}{Passed to \code{jk1weights}, \code{jknweights},
    \code{brrweights}, \code{bootweights}, \code{subbootweights}, or \code{mrbweights}.}
  \item{compress}{Use a compressed representation of the replicate
    weights matrix.}
  \item{mse}{if \code{TRUE}, compute variances from sums of squares around
    the point estimate, rather than the mean of the replicates}
}

\value{
  Object of class \code{svyrep.design}.
}
\references{
  Canty AJ, Davison AC. (1999) Resampling-based variance
  estimation for labour force surveys. The Statistician 48:379-391

  Judkins, D. (1990), "Fay's Method for Variance Estimation," Journal of Official Statistics, 6, 223-239.

  Preston J. (2009) Rescaled bootstrap for stratified multistage sampling. Survey Methodology 35(2) 227-234

  Rao JNK, Wu CFJ. Bootstrap inference for sample surveys. Proc Section
  on Survey Research Methodology. 1993 (866--871)
}

\seealso{\code{\link{brrweights}}, \code{\link{svydesign}},
  \code{\link{svrepdesign}}, \code{\link{bootweights}}, \code{\link{subbootweights}}, \code{\link{mrbweights}}
}

\examples{
data(scd)
scddes<-svydesign(data=scd, prob=~1, id=~ambulance, strata=~ESA,
nest=TRUE, fpc=rep(5,6))
scdnofpc<-svydesign(data=scd, prob=~1, id=~ambulance, strata=~ESA,
nest=TRUE)

# convert to BRR replicate weights
scd2brr <- as.svrepdesign(scdnofpc, type="BRR")
scd2fay <- as.svrepdesign(scdnofpc, type="Fay",fay.rho=0.3)
# convert to JKn weights 
scd2jkn <- as.svrepdesign(scdnofpc, type="JKn")

# convert to JKn weights with finite population correction
scd2jknf <- as.svrepdesign(scddes, type="JKn")

## with user-supplied hadamard matrix
scd2brr1 <- as.svrepdesign(scdnofpc, type="BRR", hadamard.matrix=paley(11))

svyratio(~alive, ~arrests, design=scd2brr)
svyratio(~alive, ~arrests, design=scd2brr1)
svyratio(~alive, ~arrests, design=scd2fay)
svyratio(~alive, ~arrests, design=scd2jkn)
svyratio(~alive, ~arrests, design=scd2jknf)

data(api)
## one-stage cluster sample
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
## convert to JK1 jackknife
rclus1<-as.svrepdesign(dclus1)
## convert to bootstrap
bclus1<-as.svrepdesign(dclus1,type="bootstrap", replicates=100)

svymean(~api00, dclus1)
svytotal(~enroll, dclus1)

svymean(~api00, rclus1)
svytotal(~enroll, rclus1)

svymean(~api00, bclus1)
svytotal(~enroll, bclus1)

dclus2<-svydesign(id = ~dnum + snum, fpc = ~fpc1 + fpc2, data = apiclus2)
mrbclus2<-as.svrepdesign(dclus2, type="mrb",replicates=100)
svytotal(~api00+stype, dclus2)
svytotal(~api00+stype, mrbclus2)
}
\keyword{survey}% at least one, from doc/KEYWORDS

