\name{hadamard}
\alias{hadamard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hadamard matrices }
\description{
  Returns a Hadamard matrix of dimension larger than the argument.
  Matrices of dimension every multiple of 4 are thought to exist, but
  this function doesn't know about them, so it will return matrices
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

\seealso{\code{\link{brrweights}}}
\examples{
image(hadamard(30))

plot(sapply(2:80,function(i) ncol(hadamard(i))),type="S",
     ylab="Matrix size",xlab="n",xlim=c(1,100),ylim=c(1,100))
abline(0,1,lty=2)
lines(2:80, 2:80-(2:80 \%\% 4)+4,col="purple",type="S")
legend(c(x=10,y=80),legend=c("Actual size","Minimum possible size"),
     lty=c(1,1),col=c("black","purple"),bty="n")

}
\keyword{survey}

