\name{svydesign}
\alias{svydesign}
\alias{summary.survey.design}
\alias{print.summary.survey.design}
\alias{print.survey.design}
\alias{[.survey.design}
\alias{[<-.survey.design}
\alias{na.omit.survey.design}
\alias{na.exclude.survey.design}
\alias{na.fail.survey.design}
\alias{dim.survey.design}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Survey sample analysis.}
\description{
  Specify a complex survey design.
}
\usage{
svydesign(ids, probs, strata = NULL, variables = NULL, data = NULL,     nest = FALSE, check.strata = TRUE) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{ids}{Formula or data frame specifying cluster ids from largest level to smallest level, \code{~0} is a formula for no clusters.}
  \item{probs}{Formula or data frame specifying cluster sampling probabilities}
  \item{strata}{Formula or factor specifying strata, use \code{NULL} for no strata}
  \item{variables}{Formula or data frame specifying the variables measured in the survey. If \code{NULL}, the \code{data} argument is used.}
  \item{data}{Data frame to look up variables in the formula arguments}
  \item{nest}{If \code{TRUE}, relabel cluster ids to enforce nesting}
  \item{check.strata}{If \code{TRUE}, check that clusters are nested in strata}.
}
\details{
   When analysing data from a complex survey, observations must be weighted inversely to their sampling probabilities, and the effects of stratification and of correlation induced by cluster sampling must be incorporated in standard errors.   

   The \code{svydesign} object combines a data frame and all the survey design information needed to analyse it.  These objects are used by the survey modelling and summary functions.

The  \code{dim}, \code{"["} and \code{"[<-"} and na.action methods for \code{survey.design} objects operate on the dataframe specified by \code{variables} and ensure that the design information is properly updated to correspond to the new data frame.  With the \code{"[<-"} method the new value can be a \code{survey.design} object instead of a data frame, but only the data frame is used.

}
\value{
An object of class \code{survey.design}.
}
\references{ ~put references to the literature/web site here ~ }
\author{Thomas Lumley}


\seealso{ \code{\link{svyglm}}, \code{\link{svymean}}, \code{\link{svyvar}}, \code{\link{svytable}}, \code{\link{svyquantile}}}

\examples{
  #population
  df<-data.frame(x=rnorm(1000),z=rep(0:4,200))
  df$y<-with(df, 3+3*x*z)
  #sampling fraction
  df$p<-with(df, exp(x)/(1+exp(x)))
  #sample
  xi<-rbinom(1000,1,df$p)
  sdf<-df[xi==1,]
  
  #survey design object: independent sampling, 
  dxi<-svydesign(~0,~p,data=sdf) 
  dxi
  summary(dxi)
  svymean(sdf$x,dxi)	
  svymean(~x,dxi)
  svytable(~z, dxi)

   #cluster sampling: population
   df$id<-rep(1:250,each=4)
   df$clustp<-by(df,list(df$id),function(d) min(exp(d$x*d$z)/(1+exp(d$x*d$z))))[df$id]
   xi<-rbinom(250,1,df$clustp[4*(1:250)])
   sdf<-df[xi[df$id]==1,]
   
   #cluster sampling design
   dxi<-svydesign(~id,~clustp,data=sdf)
   
   dxi
   summary(dxi)
   svymean(~x+z,dxi)

   ## stratification
   df<-data.frame(z=rep(1:4,each=200), y=rnorm(800, rep(1:4,each=200)))
   xi<-c(sample(1:200,20), sample(201:400,20), sample(401:600,20), sample(601:800,20))
   sdf<-df[xi,]
   stratdx<-svydesign(id=~0,prob=~0,strata=~z,data=sdf)
   unstrat<-svydesign(id=~0,prob=~0,data=sdf)
   stratdx
   unstrat
   summary(stratdx)

   svymean(~y, stratdx)  ##higher precision
   svymean(~y, unstrat)

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
\keyword{manip}% __ONLY ONE__ keyword per line
