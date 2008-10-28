\name{pchisqsum}
\alias{pchisqsum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distribution of quadratic forms }
\description{
  The distribution of a quadratic form in p standard Normal variables is
  a linear combination of p chi-squared distributions with 1df.
}
\usage{
pchisqsum(x, df, a, lower.tail = TRUE, method = c("satterthwaite", "integration","saddlepoint"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Observed values}
  \item{df}{Vector of degrees of freedom}
  \item{a}{Vector of coefficients }
  \item{lower.tail}{ lower or upper  tail? }
  \item{method}{See Details below}
}

\value{
Vector of cumulative probabilities
}

\details{
  The \code{"satterthwaite"} method uses Satterthwaite's approximation,
  and this is also used as a fallback for the other methods.

  \code{"integration"} inverts the characteristic function
  numerically. This is relatively slow, and not reliable for p-values
  below about 1e-5 in the upper tail, but is highly accurate for moderate p-values

  \code{"saddlepoint"} uses a saddlepoint approximation when
  \code{x>1.2*sum(lambda)} and the Satterthwaite approximation for
  smaller \code{x}. This is fast and is accurate in the upper tail, where
  accuracy is important.

  }

\references{
Davies RB (1973). "Numerical inversion of a characteristic function"
Biometrika 60:415-7

Kuonen D (1999) Saddlepoint Approximations for Distributions of Quadratic Forms in Normal Variables
 Biometrika, Vol. 86, No. 4 (Dec., 1999), pp. 929-935 

}
\seealso{\code{\link{pchisq}}}
\examples{
x <- 5*rnorm(1001)^2+rnorm(1001)^2
x.thin<-sort(x[1+(0:100)*10])
p.invert<-pchisqsum(x.thin,df=c(1,1),a=c(5,1),method="int" ,lower=FALSE)
p.satt<-pchisqsum(x.thin,df=c(1,1),a=c(5,1),method="satt",lower=FALSE)
p.sadd<-pchisqsum(x.thin,df=c(1,1),a=c(5,1),method="sad",lower=FALSE)

plot(p.invert, p.satt,type="l",log="xy")
abline(0,1,lty=2,col="purple")
plot(p.invert, p.sadd,type="l",log="xy")
abline(0,1,lty=2,col="purple")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{distribution}% __ONLY ONE__ keyword per line
