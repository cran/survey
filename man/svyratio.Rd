\name{svyratio}
\alias{svyratio}
\alias{svrepratio}
\alias{print.svyratio}
\alias{predict.svyratio}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Ratio estimation}
\description{
Ratio estimation and estimates of totals based on ratios for complex
survey samples.
}
\usage{
svyratio(numerator, denominator, design)
svrepratio(numerator, denominator, design)
\method{predict}{svyratio}(object, total, se=TRUE,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{numerator}{formula, expression, or data frame giving numerator variable(s)}
  \item{denominator}{formula, expression, or data frame giving
    denominator variable(s) }
  \item{design}{from \code{svydesign} for \code{svyratio}, from
    \code{svrepdesign} for \code{svrepratio}}
  \item{object}{result of \code{svyratio}}
  \item{total}{vector of population totals for the denominator variables in
    \code{object}}
  \item{se}{Return standard errors?}
  \item{...}{Other unused arguments to \code{predict}}
}

\value{
\code{svyratio} returns an object of class \code{svyratio}. The
\code{predict} method returns a matrix of population totals and
optionally a matrix of standard errors.
}
\references{Levy and Lemeshow. "Sampling of Populations" (3rd edition). Wiley}
\author{Thomas Lumley}

\seealso{\code{\link{svydesign}} }

\examples{
data(scd)

## survey design objects
scddes<-svydesign(data=scd, prob=~1, id=~ambulance, strata=~ESA,
nest=TRUE, fpc=rep(5,6))
scdnofpc<-svydesign(data=scd, prob=~1, id=~ambulance, strata=~ESA,
nest=TRUE)

# convert to BRR replicate weights
scd2brr <- as.svrepdesign(scdnofpc, type="BRR")

# use BRR replicate weights from Levy and Lemeshow
repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
c(0,1,0,1,1,0))
scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)

# ratio estimates
svyratio(~alive, ~arrests, design=scddes)
svyratio(~alive, ~arrests, design=scdnofpc)
svrepratio(~alive, ~arrests, design=scd2brr)
svrepratio(~alive, ~arrests, design=scdrep)

}
\keyword{survey}% at least one, from doc/KEYWORDS

