\name{svysmooth}
\alias{svysmooth}
\alias{svysmooth.default}
\alias{plot.svysmooth}
\alias{print.svysmooth}
\alias{lines.svysmooth}
\alias{make.panel.svysmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scatterplot smoothing and density estimation}
\description{
Scatterplot smoothing and density estimation for probability-weighted
data. 
}
\usage{
svysmooth(formula, design, ...)
\method{svysmooth}{default}(formula, design, method = c("locpoly", "quantreg"), 
    bandwidth = NULL, quantile, df = 4, ...)
\method{plot}{svysmooth}(x, which=NULL, type="l", xlabs=NULL, ylab=NULL,...)
\method{lines}{svysmooth}(x,which=NULL,...)
make.panel.svysmooth(design,bandwidth=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{One-sided formula for density estimation, two-sided for smoothing}
  \item{design}{Survey design object}
  \item{method}{local polynomial smoothing for the mean or regression
    splines for quantiles}
  \item{bandwidth}{Smoothing bandwidth for "locpoly" or \code{NULL} for automatic choice}
  \item{quantile}{quantile to be estimated for "quantreg"}
  \item{df}{Degrees of freedom for "quantreg"}
  \item{which}{Which plots to show (default is all)}
  \item{type}{as for \code{plot}}
  \item{xlabs}{Optional vector of x-axis labels}
  \item{ylab}{Optional y-axis label}
  \item{\dots}{More arguments}
  \item{x}{Object of class \code{svysmooth}}
}
\details{
\code{svysmooth} does one-dimensional smoothing. If \code{formula} has
multiple predictor variables a separate one-dimensional smooth is
performed for each one. 
  
For \code{method="locpoly"} the extra arguments are passed to
\code{locpoly} from the KernSmooth package, for
\code{method="quantreg"} they are passed to \code{rq} from the
quantreg package.  The automatic choice of bandwidth for
\code{method="locpoly"} uses the default settings for \code{dpik} and
\code{dpill} in the KernSmooth package.

\code{make.panel.svysmooth()} makes a function that plots points and
draws a weighted smooth curve through them, a weighted replacement for
\code{\link{panel.smooth}} that can be passed to functions such as
\code{\link{termplot}} or \code{\link{plot.lm}}.  The resulting function has a \code{span} argument that will set the bandwidth; if this is not specified the automatic choice will be used.
}
\value{
An object of class \code{svysmooth}, a list of lists, each with \code{x} and \code{y} components.
}

\seealso{\code{\link{svyhist}} for histograms}
\examples{
 data(api)
 dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)

 smth<-svysmooth(api00~api99+ell,dstrat)
 dens<-svysmooth(~api99, dstrat,bandwidth=30)
 dens1<-svysmooth(~api99, dstrat)
 qsmth<-svysmooth(api00~ell,dstrat, quantile=0.75, df=3,method="quantreg")

 plot(smth)
 plot(smth, which="ell",lty=2,ylim=c(500,900))
 lines(qsmth, col="red")

 svyhist(~api99,design=dstrat)
 lines(dens,col="purple",lwd=3)
 lines(dens1, col="forestgreen",lwd=2)

 m<-svyglm(api00~sin(api99/100)+stype, design=dstrat)
 termplot(m, data=model.frame(dstrat), partial.resid=TRUE, se=TRUE,
  smooth=make.panel.svysmooth(dstrat))
}
\keyword{hplot}% at least one, from doc/KEYWORDS
\keyword{survey}

