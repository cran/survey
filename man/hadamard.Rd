\name{hadamard}
\alias{hadamard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hadamard matrices }
\description{
  Returns a Hadamard matrix of dimension larger than the argument.
  Matrices of dimension every multiple of 4 are thought to exist, but
  this function doesn't know about all of them, so it will sometimes return matrices
  that are larger than necessary.
}
\usage{
hadamard(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{lower bound for size }
}

\value{
A Hadamard matrix
}

\note{Strictly speaking, a Hadamard matrix has entries +1 and -1 rather
  than 1 and 0, so \code{2*hadamard(n)-1} is a Hadamard matrix}
\references{
  Sloane NJA. A Library of Hadamard Matrices \url{http://www.research.att.com/~njas/hadamard/}
  }
\seealso{\code{\link{brrweights}}}
\examples{

par(mfrow=c(2,2))
## Sylvester-type
image(hadamard(63),main="Sylvester: 64")
## Paley-type
image(hadamard(57),main="Paley: 60")
## hybrid
image(hadamard(23),main="Stored: 24")
image(hadamard(90),main="Constructed: 96=4x24")

par(mfrow=c(1,1))
plot(2:150,sapply(2:150,function(i) ncol(hadamard(i))),type="S",
     ylab="Matrix size",xlab="n",xlim=c(1,150),ylim=c(1,150))
abline(0,1,lty=3)
lines(2:150, 2:150-(2:150 \%\% 4)+4,col="purple",type="S",lty=2)
legend(c(x=10,y=140),legend=c("Actual size","Minimum possible size"),
     col=c("black","purple"),bty="n",lty=c(1,2))

}
\keyword{survey}

