\name{ftable.svystat}
\alias{ftable.svystat}
\alias{ftable.svrepstat}
\alias{ftable.svyby}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lay out tables of survey statistics}
\description{
Reformat the output of survey computations to a table.  
}
\usage{
\method{ftable}{svystat}(x, rownames,...)
\method{ftable}{svrepstat}(x, rownames,...)
\method{ftable}{svyby}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Output of functions such as \code{svymean},\code{svrepmean}, \code{svyby}}
  \item{rownames}{List of vectors of strings giving dimension names for
    the resulting table (see examples)}
  \item{...}{Arguments for future expansion}
}
\value{
 An object of class \code{"ftable"}
}

\seealso{ \code{\link{ftable}}}
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

a<-svymean(~interaction(stype,comp.imp), design=dclus1)
b<-ftable(a, rownames=list(stype=c("E","H","M"),comp.imp=c("No","Yes")))
b

a<-svymean(~interaction(stype,comp.imp), design=dclus1, deff=TRUE)
b<-ftable(a, rownames=list(stype=c("E","H","M"),comp.imp=c("No","Yes")))
round(100*b,1)

rclus1<-as.svrepdesign(dclus1)
a<-svytotal(~interaction(stype,comp.imp), design=rclus1)
b<-ftable(a, rownames=list(stype=c("E","H","M"),comp.imp=c("No","Yes")))
b
round(b)

a<-svyby(~api99 + api00, ~stype + sch.wide, rclus1, svymean, keep.var=TRUE)
ftable(a)
print(ftable(a),digits=2)

b<-svyby(~api99 + api00, ~stype + sch.wide, rclus1, svymean, keep.var=TRUE, deff=TRUE)
print(ftable(b),digits=2)

d<-svyby(~api99 + api00, ~stype + sch.wide, rclus1, svymean, keep.var=TRUE, vartype=c("se","cvpct"))
round(ftable(d),1)

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
