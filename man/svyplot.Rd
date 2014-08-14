\name{svyplot}
\alias{svyplot}
\alias{svyplot.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots for survey data }
\description{
  Because observations in survey samples may represent very different
  numbers of units in the population ordinary plots can be misleading.
  The \code{svyplot} function produces scatterplots adjusted in various ways
  for sampling weights. 
}
\usage{
svyplot(formula, design,...)
\method{svyplot}{default}(formula, design, style = c("bubble", "hex", "grayhex","subsample","transparent"),
sample.size = 500, subset = NULL, legend = 1, inches = 0.05,
amount=NULL, basecol="black",
alpha=c(0, 0.8),xbins=30,...)
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
    jittering to use in subsample plots, or \code{NULL} for the default
    amount}
  \item{basecol}{base color for transparent plots, or a function to
    compute the color (see below), or color for bubble plots}
  \item{alpha}{minimum and maximum opacity for transparent plots}
\item{xbins}{Number of (x-axis) bins for hexagonal binning}
  \item{\dots}{Passed to \code{plot} methods}
}
\details{
  Bubble plots are scatterplots with circles whose area is proportional
  to the sampling weight.  The two "hex" styles produce hexagonal
  binning scatterplots, and require the \code{hexbin} package from
  Bioconductor. The "transparent" style plots points with opacity
  proportional to sampling weight.

  
  The \code{subsample} method uses the sampling weights to create a
  sample from approximately the population distribution and passes this to \code{\link{plot}}
  
  Bubble plots are suited to small surveys, hexagonal binning and
  transparency to large surveys where plotting all the points would
  result in too much overlap.

  \code{basecol} can be a function taking one data frame argument, which
  will be passed the data frame of variables from the survey object.
  This could be memory-intensive for large data sets.
}
\value{
  None
}

\seealso{
  \code{\link{symbols}} for other options (such as colour) for bubble
  plots.

  \code{\link{svytable}} for plots of discrete data.
}
\examples{
data(api)
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)

svyplot(api00~api99, design=dstrat, style="bubble")
svyplot(api00~api99, design=dstrat, style="transparent",pch=19)

## these two require the hexbin package 
svyplot(api00~api99, design=dstrat, style="hex", xlab="1999 API",ylab="2000 API")
svyplot(api00~api99, design=dstrat, style="grayhex",legend=0)


dclus2<-svydesign(id=~dnum+snum,  weights=~pw,
                    data=apiclus2, fpc=~fpc1+fpc2)
svyplot(api00~api99, design=dclus2, style="subsample")
svyplot(api00~api99, design=dclus2, style="subsample",
          amount=list(x=25,y=25))

svyplot(api00~api99, design=dstrat,
  basecol=function(df){c("goldenrod","tomato","sienna")[as.numeric(df$stype)]},
  style="transparent",pch=19,alpha=c(0,1))
legend("topleft",col=c("goldenrod","tomato","sienna"), pch=19, legend=c("E","H","M"))

## For discrete data, estimate a population table and plot the table.
plot(svytable(~sch.wide+comp.imp+stype,design=dstrat))
fourfoldplot(svytable(~sch.wide+comp.imp+stype,design=dstrat,round=TRUE))


## To draw on a hexbin plot you need grid graphics, eg,
library(grid)
h<-svyplot(api00~api99, design=dstrat, style="hex", xlab="1999 API",ylab="2000 API")
s<-svysmooth(api00~api99,design=dstrat)
grid.polyline(s$api99$x,s$api99$y,vp=h$plot.vp@hexVp.on,default.units="native", 
   gp=gpar(col="red",lwd=2))
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line
