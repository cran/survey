\name{brrweights}
\alias{jk1weights}
\alias{jknweights}
\alias{brrweights}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Compute replicate weights }
\description{
  Compute replicate weights from a survey design. These functions are
  usually called from \code{\link{as.svrepdesign}} rather than directly
  by the user.
}
\usage{
brrweights(strata, psu, match = NULL,
              small = c("fail","split","merge"),
              large = c("split", "merge", "fail"),
              fay.rho=0, only.weights=FALSE,
              compress=TRUE, hadamard.matrix=NULL)
jk1weights(psu,fpc=NULL,
              fpctype=c("population","fraction","correction"),
              compress=TRUE)
jknweights(strata,psu, fpc=NULL,
              fpctype=c("population","fraction","correction"),
              compress=TRUE,
              lonely.psu=getOption("survey.lonely.psu"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{strata}{Stratum identifiers }
  \item{psu}{PSU (cluster) identifier }
  \item{match}{Optional variable to use in matching. }
  \item{small}{How to handle strata with only one PSU}
  \item{large}{How to handle strata with more than two PSUs}
  \item{fpc}{Optional population (stratum) size or finite population correction }
  \item{fpctype}{How \code{fpc} is coded.}
  \item{fay.rho}{Parameter for Fay's extended BRR method}
  \item{only.weights}{If \code{TRUE} return only the matrix of
    replicate weights}
  \item{compress}{If \code{TRUE}, store the replicate weights in
    compressed form}
  \item{hadamard.matrix}{Optional user-supplied Hadamard matrix for
    \code{brrweights}}
  \item{lonely.psu}{Handling of non-certainty single-PSU strata}
}
\details{
  JK1 and JKn are  jackknife schemes for unstratified and stratified
  designs respectively.  The finite population correction may be
  specified as a single number, a vector with one entry per stratum, or
  a vector with one entry per observation (constant within strata). 
  When \code{fpc} is a vector with one entry per stratum it may not have
  names that differ from the stratum identifiers (it may have no names,
  in which case it must be in the same order as
  \code{unique(strata)}). To specify population stratum sizes use
  \code{fpctype="population"}, to specify sampling fractions use
  \code{fpctype="fraction"} and to specify the correction directly use
  \code{fpctype="correction"}
  
  The only reason not to use \code{compress=TRUE} is that it is new and
  there is a greater possibility of bugs.  It reduces the number of
  rows of the replicate weights matrix from the number of observations
  to the number of PSUs.
  
  In BRR variance estimation each stratum is split in two to give
  half-samples. Balanced replicated weights are needed, where
  observations in two different strata end up in the same half stratum
  as often as in different half-strata.BRR, strictly speaking, is
  defined only when each stratum has exactly
  two PSUs.  A stratum with one PSU can be merged with another such
  stratum, or can be split to appear in both half samples with half
  weight.  The latter approach is appropriate for a PSU that was
  deterministically sampled.
  
  A stratum with more than two PSUs can be split into multiple smaller
  strata each with two PSUs or the PSUs can be merged to give two
  superclusters within the stratum.
  
  When merging small strata or grouping PSUs in large strata the
  \code{match} variable is used to sort PSUs before merging, to give
  approximate matching on this variable.
  
  If you want more control than this you should probably construct your
  own weights using the Hadamard matrices produced by \code{\link{hadamard}}
  

}
\value{
For \code{brrweights} with \code{only.weights=FALSE} a list with elements
  \item{weights}{two-column matrix indicating the weight for each
    half-stratum in one particular set  of split samples}
  \item{wstrata}{New stratum variable incorporating merged or split strata}
  \item{strata}{Original strata for distinct PSUs}
  \item{psu}{Distinct PSUs}
  \item{npairs}{Dimension of Hadamard matrix used in BRR construction}
  \item{sampler}{function returning replicate weights}
  \item{compress}{Indicates whether the \code{sampler} returns per PSU
    or per observation weights}
  For \code{jk1weights} and \code{jknweights} a data frame of replicate
  weights and the \code{scale} and \code{rscale} arguments to \code{\link{svrVar}}.
}
\references{Levy and Lemeshow "Sampling of Populations". Wiley.
  
  Shao and Tu "The Jackknife and Bootstrap". Springer.
}

\seealso{\code{\link{hadamard}}, \code{\link{as.svrepdesign}},
  \code{\link{svrVar}}, \code{\link{surveyoptions}}}

\examples{
data(scd)
scdnofpc<-svydesign(data=scd, prob=~1, id=~ambulance, strata=~ESA,
nest=TRUE)

## convert to BRR replicate weights
scd2brr <- as.svrepdesign(scdnofpc, type="BRR")
svymean(~alive, scd2brr)
svyratio(~alive, ~arrests, scd2brr)

## with user-supplied hadamard matrix
scd2brr1 <- as.svrepdesign(scdnofpc, type="BRR", hadamard.matrix=paley(11))
svymean(~alive, scd2brr1)
svyratio(~alive, ~arrests, scd2brr1)


}
\keyword{survey}% at least one, from doc/KEYWORDS

