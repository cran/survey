\name{barplot.svystat}
\alias{barplot.svystat}
\alias{barplot.svrepstat}
\alias{barplot.svyby}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Barplots }
\description{
Draws a barplot based on results from a survey analysis. The default
barplot method already works for results from \code{\link{svytable}}.
}
\usage{
\method{barplot}{svystat}(height, ...)

\method{barplot}{svrepstat}(height, ...)

\method{barplot}{svyby}(height,beside=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{height}{Analysis result }
  \item{beside}{Grouped, rather than stacked, bars}
  \item{\dots}{ Arguments to \code{\link{barplot}} }
}


\examples{

data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

a<-svymean(~stype, dclus1)
barplot(a)
barplot(a, names.arg=c("Elementary","High","Middle"), col="purple", main="Proportions of school level")

b<-svyby(~enroll+api.stu, ~stype, dclus1, svymean)
barplot(b,beside=TRUE,legend=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{hplot}% __ONLY ONE__ keyword per line
