\name{svycontrast}
\alias{svycontrast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Linear and nonlinearconstrasts of survey statistics }
\description{
Computes linear or nonlinear contrasts of estimates produced by survey
functions (or any object with \code{coef} and \code{vcov} methods).
}
\usage{
svycontrast(stat, contrasts, add=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stat}{object of class \code{svrepstat} or \code{svystat} }
  \item{contrasts}{A vector or list of vectors of coefficients, or a
    call or list of calls }
  \item{add}{keep all the coefficients of the input in the output?}
  \item{\dots}{For future expansion}
}
\value{
Object of class \code{svrepstat} or \code{svystat} or \code{svyvar}
}

\details{
If \code{contrasts} is a list, the element names are used as
names for the returned statistics.
  
If an element of \code{contrasts} is shorter than \code{coef(stat)} and has names, the
names are used to match up the vectors and the remaining elements of
\code{contrasts} are assumed to be zero. If the names are not legal
variable names (eg \code{0.1}) they must be quoted (eg \code{"0.1"})

If \code{contrasts} is a \code{"call"} or list of \code{"call"s}, and
\code{stat} is a \code{svrepstat} object including replicates, the
replicates are transformed and used to compute the variance. If
\code{stat} is a \code{svystat} object or a \code{svrepstat} object
without replicates, the delta-method is used to compute variances, and
the calls must use only functions that \code{\link{deriv}} knows how
to differentiate. If the names are not legal variable names they must
be quoted with backticks (eg \code{`0.1`}).

If \code{stats} is a \code{svyvar} object, the estimates are elements of a matrix and the names are the row and column names pasted together with a colon. 
  }
\seealso{\code{\link{regTermTest}}, \code{svyglm}}
\examples{
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

a <- svytotal(~api00+enroll+api99, dclus1)
svycontrast(a, list(avg=c(0.5,0,0.5), diff=c(1,0,-1)))
## if contrast vectors have names, zeroes may be omitted
svycontrast(a, list(avg=c(api00=0.5,api99=0.5), diff=c(api00=1,api99=-1)))

## nonlinear contrasts
svycontrast(a, quote(api00/api99))
svyratio(~api00, ~api99, dclus1)

## Example: standardised skewness coefficient
moments<-svymean(~I(api00^3)+I(api00^2)+I(api00), dclus1)
svycontrast(moments, 
quote((`I(api00^3)`-3*`I(api00^2)`*`I(api00)`+ 3*`I(api00)`*`I(api00)`^2-`I(api00)`^3)/
      (`I(api00^2)`-`I(api00)`^2)^1.5))

## Example: geometric means
## using delta method
meanlogs <- svymean(~log(api00)+log(api99), dclus1)
svycontrast(meanlogs,
    list(api00=quote(exp(`log(api00)`)), api99=quote(exp(`log(api99)`))))

## using delta method
rclus1<-as.svrepdesign(dclus1)
meanlogs <- svymean(~log(api00)+log(api99), rclus1)
svycontrast(meanlogs,
    list(api00=quote(exp(`log(api00)`)),
api99=quote(exp(`log(api99)`))))

## why is add=TRUE useful?

(totals<-svyby(~enroll,~stype,design=dclus1,svytotal,covmat=TRUE))
totals1<-svycontrast(totals, list(total=c(1,1,1)), add=TRUE)

svycontrast(totals1, list(quote(E/total), quote(H/total), quote(M/total)))


totals2<-svycontrast(totals, list(total=quote(E+H+M)), add=TRUE)
all.equal(as.matrix(totals1),as.matrix(totals2))

## more complicated svyby
means <- svyby(~api00+api99, ~stype+sch.wide, design=dclus1, svymean,covmat=TRUE)
svycontrast(means, quote(`E.No:api00`-`E.No:api99`))
svycontrast(means, quote(`E.No:api00`/`E.No:api99`))

## transforming replicates
meanlogs_r <- svymean(~log(api00)+log(api99), rclus1, return.replicates=TRUE)
svycontrast(meanlogs_r,
    list(api00=quote(exp(`log(api00)`)), api99=quote(exp(`log(api99)`))))

## converting covariances to correlations
vmat <-svyvar(~api00+ell,dclus1)
print(vmat,cov=TRUE)
cov2cor(as.matrix(vmat))[1,2]
svycontrast(vmat, quote(`api00:ell`/sqrt(`api00:api00`*`ell:ell`)))

}
\keyword{survey}% at least one, from doc/KEYWORDS

