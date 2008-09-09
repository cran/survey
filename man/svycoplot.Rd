\name{svycoplot}
\alias{svycoplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conditioning plots of survey data }
\description{
Draws conditioned scatterplots ('Trellis' plots) of survey data using
hexagonal binning or transparency.
}
\usage{
svycoplot(formula, design, style = c("hexbin", "transparent"), basecol =
"black", alpha = c(0, 0.8),hexscale=c("relative","absolute"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A graph formula suitable for \code{\link{xyplot}}}
  \item{design}{A survey design object }
  \item{style}{Hexagonal binning or transparent color?}
  \item{basecol}{The fully opaque 'base' color for creating transparent
    colors. This may also be a function; see \code{\link{svyplot}} for details}
  \item{alpha}{Minimum and maximum opacity }
  \item{hexscale}{Scale hexagons separate for each panel (relative) or
    across all panels (absolute)}
  \item{\dots}{Other arguments passed to \code{grid.hexagons} or \code{\link{xyplot}} }
}

\value{
  An object of class \code{trellis}
}
\note{
As with all 'Trellis' graphs, this function creates an object but does
  not draw the graph. When used inside a function or non-interactively
  you need to \code{print()} the result to create the graph.
}
\seealso{\code{\link{svyplot}}}
\examples{
data(api)
dclus2<-svydesign(id=~dnum+snum,  weights=~pw,
                    data=apiclus2, fpc=~fpc1+fpc2)

svycoplot(api00~api99|sch.wide*comp.imp, design=dclus2, style="hexbin")
svycoplot(api00~api99|sch.wide*comp.imp, design=dclus2, style="hexbin", hexscale="absolute")

svycoplot(api00~api99|sch.wide, design=dclus2, style="trans")

svycoplot(api00~meals|stype,design=dclus2,
  style="transparent",
  basecol=function(d) c("darkred","purple","forestgreen")[as.numeric(d$stype)],
  alpha=c(0,1)) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{hplot}% __ONLY ONE__ keyword per line
