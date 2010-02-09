\name{svyprcomp}
\Rdversion{1.1}
\alias{svyprcomp}
\alias{biplot.svyprcomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampling-weighted principal component analysis
}
\description{
Computes principal components using the sampling weights. 
}
\usage{
svyprcomp(formula, design, center = TRUE, scale. = FALSE, tol = NULL, scores = FALSE, ...)
\method{biplot}{svyprcomp}(x, cols=c("black","darkred"),xlabs=NULL,weight=c("transparent","scaled","none"),
                           max.alpha=0.5,max.cex=0.5,xlim=NULL,ylim=NULL,pc.biplot=FALSE,expand=1,xlab=NULL,ylab=NULL,
                           arrow.len=0.1,
                           ...)
}
\arguments{
  \item{formula}{
model formula describing variables to be used
}
  \item{design}{
survey design object.
}
  \item{center}{
Center data before analysis?
}
  \item{scale.}{
Scale to unit variance before analysis?
}
  \item{tol}{
Tolerance for omitting components from the results; a proportion of the standard deviation of the first component.  The default is to keep all components.
}
  \item{scores}{
Return scores on each component? These are needed for \code{biplot}.
}
\item{x}{
A \code{svyprcomp} object
}
\item{cols}{
Base colors for observations and variables respectively
}
\item{xlabs}{
 Formula, or character vector, giving labels for each observation
}
\item{weight}{
 How to display the sampling weights: \code{"scaled"} changes the size of the point label, \code{"transparent"} uses opacity proportional to sampling weight, \code{"none"} changes neither. 
}
\item{max.alpha}{
 Opacity for the largest sampling weight, or for all points if \code{weight!="transparent"}
}
\item{max.cex}{
 Character size (as a multiple of \code{par("cex")}) for the largest sampling weight, or for all points if \code{weight!="scaled"}
}
\item{xlim,ylim,xlab,ylab}{Graphical parameters}
\item{expand,arrow.len}{See \code{\link{biplot}}}
\item{pc.biplot}{See \code{link{biplot.prcomp}}}
  \item{\dots}{	    
Other arguments to \code{\link{prcomp}}, or graphical parameters for \code{biplot}
}

}

\value{
\code{svyprcomp} returns an object of class \code{svyprcomp}, similar to
class \code{prcomp} but including design information

}

\seealso{
\code{\link{prcomp}}, \code{\link{biplot.prcomp}}
}
\examples{
data(api)
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)

pc <- svyprcomp(~api99+api00+ell+hsg+meals+emer, design=dclus2,scale=TRUE,scores=TRUE)
pc
biplot(pc, xlabs=~dnum, weight="none")

biplot(pc, xlabs=~dnum,max.alpha=1)

biplot(pc, weight="scaled",max.cex=1.5, xlabs=~dnum)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{hplot}% __ONLY ONE__ keyword per line
\keyword{multivariate}