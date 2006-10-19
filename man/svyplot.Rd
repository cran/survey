\name{svyplot}
\alias{svyplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots for survey data }
\description{
  Because observations in survey samples may represent very different
  numbers of units in the population ordinary plots can be misleading.
  The \code{svyplot} function produces plots adjusted in various ways
  for sampling weights.
}
\usage{
svyplot(formula, design, style = c("bubble", "hex", "grayhex","subsample"),
sample.size = 500, subset = NULL, legend = 1, inches = 0.05, amount=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A model formula}
  \item{design}{ A survey object (svydesign or svrepdesign)}
  \item{style}{See Details below}
  \item{sample.size}{For \code{style="subsample"}}
  \item{subset}{expression using variables in the design object}
  \item{legend}{For \code{style="hex"} or \code{"grayhex"}}
  \item{inches}{Scale for bubble plots}
  \item{amount}{list with \code{x} and \code{y} components for amount of
    jittering to use in subsample plots, or \code{NULL} for the default amount}
  \item{\dots}{Passed to \code{plot} methods}
}
\details{
  Bubble plots are scatterplots with circles whose area is proportional
  to the sampling weight.  The two "hex" styles produce hexagonal
  binning scatterplots, and require the \code{hexbin} package from
  Bioconductor.
  
  The \code{subsample} method uses the sampling weights to create a
  sample from approximately the population distribution and passes this to \code{\link{plot}}
  
  Bubble plots are suited to small surveys, hexagonal binning and
  subsampling to large surveys where plotting all the points would
  result in too much overlap.
}
\value{
  None
}

\seealso{
  \code{\link{symbols}} for other options (such as colour) for bubble plots.
}
\examples{
data(api)
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)

svyplot(api00~api99, design=dstrat, style="bubble")
\dontrun{
## these two require the hexbin package from Bioconductor
svyplot(api00~api99, design=dstrat, style="hex", xlab="1999 API",ylab="2000 API")
svyplot(api00~api99, design=dstrat, style="grayhex",legend=0)
}

dclus2<-svydesign(id=~dnum+snum,  weights=~pw,
                    data=apiclus2, fpc=~fpc1+fpc2)
svyplot(api00~api99, design=dclus2, style="subsample")
svyplot(api00~api99, design=dclus2, style="subsample",
          amount=list(x=25,y=25))

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line
