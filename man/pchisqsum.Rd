\name{pchisqsum}
\alias{pchisqsum}
\alias{pFsum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distribution of quadratic forms }
\description{
  The distribution of a quadratic form in p standard Normal variables is
  a linear combination of p chi-squared distributions with 1df.   When
  there is uncertainty about the variance, a reasonable model for the
  distribution is a linear combination of F distributions with the same
  denominator. 

}
\usage{
pchisqsum(x, df, a, lower.tail = TRUE, 
   method = c("satterthwaite", "integration","saddlepoint"))
pFsum(x, df, a, ddf=Inf,lower.tail = TRUE, 
   method = c("saddlepoint","integration","satterthwaite"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Observed values}
  \item{df}{Vector of degrees of freedom}
  \item{a}{Vector of coefficients }
  \item{ddf}{Denominator degrees of freedom}
  \item{lower.tail}{ lower or upper  tail? }
  \item{method}{See Details below}
  \item{\dots}{arguments to \code{pchisqsum}}
}

\value{
Vector of cumulative probabilities
}

\details{
  The \code{"satterthwaite"} method uses Satterthwaite's
  approximation, and this is also used as a fallback for the other
  methods. The accuracy is usually good, but is more variable depending
  on \code{a} than the other methods and is anticonservative in the
  extreme tail. The Satterthwaite approximation requires all \code{a>0}.

  \code{"integration"} requires the \code{CompQuadForm} package. For
  \code{pchisqsum} it uses Farebrother's algorithm if all
  \code{a>0}. For \code{pFsum} or when some \code{a<0} it inverts the
  characteristic function using the algorithm of Davies (1980).  If
  the \code{CompQuadForm} package is not present, a warning is given
  and the saddlepoint approximation is used.  These algorithms are not
  accurate for very large \code{x} or when some \code{a} are close to
  zero: a warning is given if the relative error bound is more than
  10\% of the result.
  
  \code{"saddlepoint"} uses Kuonen's saddlepoint approximation. This
  is accurate even very far out in the upper tail or with some
  \code{a=0} and does not require any additional packages. It is
  implemented in pure R and so is slower than the \code{"integration"}
  method.

  The distribution in \code{pFsum} is standardised so that a likelihood
  ratio test can use the same \code{x} value as in \code{pchisqsum}.
  That is, the linear combination of chi-squareds is multiplied by
  \code{ddf} and then divided by an independent chi-squared with
  \code{ddf} degrees of freedom.
  
  }

\references{
Davies RB (1973). "Numerical inversion of a characteristic function"
Biometrika 60:415-7

Davies RB (1980) "Algorithm AS 155: The Distribution of a Linear Combination of chi-squared Random Variables"
Applied Statistics,Vol. 29, No. 3 (1980), pp. 323-333

P. Duchesne, P. Lafaye de Micheaux (2010) "Computing the distribution
of quadratic forms: Further comparisons between the Liu-Tang-Zhang
approximation and exact methods", Computational Statistics and Data
Analysis, Volume 54, (2010), 858-862

Farebrother R.W. (1984)  "Algorithm AS 204: The distribution of a
Positive Linear Combination of chi-squared random variables". Applied
Statistics Vol. 33, No. 3 (1984), p. 332-339

Kuonen D (1999) Saddlepoint Approximations for Distributions of
Quadratic Forms in Normal Variables. Biometrika, Vol. 86, No. 4
(Dec., 1999), pp. 929-935 

}
\seealso{\code{\link{pchisq}}}
\examples{
x <- 2.7*rnorm(1001)^2+rnorm(1001)^2+0.3*rnorm(1001)^2
x.thin<-sort(x)[1+(0:50)*20]
p.invert<-pchisqsum(x.thin,df=c(1,1,1),a=c(2.7,1,.3),method="int" ,lower=FALSE)
p.satt<-pchisqsum(x.thin,df=c(1,1,1),a=c(2.7,1,.3),method="satt",lower=FALSE)
p.sadd<-pchisqsum(x.thin,df=c(1,1,1),a=c(2.7,1,.3),method="sad",lower=FALSE)

plot(p.invert, p.satt,type="l",log="xy")
abline(0,1,lty=2,col="purple")
plot(p.invert, p.sadd,type="l",log="xy")
abline(0,1,lty=2,col="purple")

pchisqsum(20, df=c(1,1,1),a=c(2.7,1,.3), lower.tail=FALSE,method="sad")
pFsum(20, df=c(1,1,1),a=c(2.7,1,.3), ddf=49,lower.tail=FALSE,method="sad")
pFsum(20, df=c(1,1,1),a=c(2.7,1,.3), ddf=1000,lower.tail=FALSE,method="sad")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{distribution}% __ONLY ONE__ keyword per line
