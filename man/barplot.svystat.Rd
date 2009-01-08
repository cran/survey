\name{barplot.svystat}
\alias{barplot.svystat}
\alias{barplot.svrepstat}
\alias{barplot.svyby}
\alias{dotchart}
\alias{dotchart.svystat}
\alias{dotchart.svrepstat}
\alias{dotchart.svyby}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Barplots and Dotplots }
\description{
Draws a barplot or dotplot based on results from a survey analysis. The default
barplot method already works for results from \code{\link{svytable}}.
}
\usage{
\method{barplot}{svystat}(height, ...)
\method{barplot}{svrepstat}(height, ...)
\method{barplot}{svyby}(height,beside=TRUE, ...)

\method{dotchart}{svystat}(x,...,pch=19)
\method{dotchart}{svrepstat}(x,...,pch=19)
\method{dotchart}{svyby}(x,...,pch=19)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{height,x}{Analysis result }
  \item{beside}{Grouped, rather than stacked, bars}
  \item{\dots}{ Arguments to \code{\link{barplot}} or \code{dotchart} }
  \item{pch}{Overrides the default in \code{dotchart.default}}
}


\examples{

data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

a<-svymean(~stype, dclus1)
barplot(a)
barplot(a, names.arg=c("Elementary","High","Middle"), col="purple", main="Proportions of school level")

b<-svyby(~enroll+api.stu, ~stype, dclus1, svymean)
barplot(b,beside=TRUE,legend=TRUE)
dotchart(b)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{hplot}% __ONLY ONE__ keyword per line
