\name{pchisqsum}
\alias{pchisqsum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distribution of quadratic forms }
\description{
The distribution of a quadratic form in p standard Normal variables is a linear
combination of p chi-squared distributions with 1df. This function
provides the cumulative distribution function by numerical inversion
of the characteristic function and also provides the Satterthwaite approximation.
}
\usage{
pchisqsum(x, df, a, lower.tail = TRUE,  method = c("satterthwaite", "integration"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Observed values }
  \item{df}{Vector of degrees of freedom}
  \item{a}{Vector of coefficients}
  \item{lower.tail}{lower tail or upper tail?}
  \item{method}{Satterthwaite approximation or numerical integration?}
}

\value{
 Vector of cumulative probabilities
}

\seealso{\code{\link{pchisq}}}
  
\examples{
x <- 5*rnorm(1001)^2+rnorm(1001)^2
x.thin<-sort(x[1+(0:50)*20])
p.invert<-pchisqsum(x.thin,df=c(1,1),a=c(5,1),method="integration")
p.satt<-pchisqsum(x.thin,df=c(1,1),a=c(5,1),method="satt")

plot(p.invert, p.satt,type="l")
abline(0,1,lty=2,col="purple")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

