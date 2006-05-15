\name{hadamard}
\alias{hadamard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hadamard matrices }
\description{
  Returns a Hadamard matrix of dimension larger than the argument.
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
\details{
  For most \code{n} the matrix comes from \code{\link{paley}}. The
  \eqn{36\times 36}{36x36} matrix is from Plackett and Burman (1946)
  and the \eqn{28\times 28}{28x28} is from Sloane's library of Hadamard
  matrices.
  
  Matrices of dimension every multiple of 4 are thought to exist, but
  this function doesn't know about all of them, so it will sometimes
  return matrices that are larger than necessary.  The excess is at most
  4 for \eqn{n<180}{n<180} and at most 5\% for \eqn{n>100}{n>100}.
  
}
\note{Strictly speaking, a Hadamard matrix has entries +1 and -1 rather
  than 1 and 0, so \code{2*hadamard(n)-1} is a Hadamard matrix}
\references{
  Sloane NJA. A Library of Hadamard Matrices \url{http://www.research.att.com/~njas/hadamard/}


Plackett RL, Burman JP. (1946) The Design of Optimum Multifactorial Experiments
Biometrika, Vol. 33, No. 4  pp. 305-325 
  
  Cameron PJ (2005) Hadamard Matrices
\url{http://designtheory.org/library/encyc/topics/had.pdf}. In: The
Encyclopedia of Design Theory \url{http://designtheory.org/library/encyc/}
  }
\seealso{\code{\link{brrweights}}, \code{\link{paley}}}
\examples{

par(mfrow=c(2,2))
## Sylvester-type
image(hadamard(63),main=quote("Sylvester: "*64==2^6))
## Paley-type
image(hadamard(59),main=quote("Paley: "*60==59+1))
## from NJ Sloane's library
image(hadamard(27),main=quote("Stored: "*28))
## For n=90 we get 96 rather than the minimum possible size, 92.
image(hadamard(90),main=quote("Constructed: "*96==2^3\%*\%(11+1)))

par(mfrow=c(1,1))
plot(2:150,sapply(2:150,function(i) ncol(hadamard(i))),type="S",
     ylab="Matrix size",xlab="n",xlim=c(1,150),ylim=c(1,150))
abline(0,1,lty=3)
lines(2:150, 2:150-(2:150 \%\% 4)+4,col="purple",type="S",lty=2)
legend(c(x=10,y=140),legend=c("Actual size","Minimum possible size"),
     col=c("black","purple"),bty="n",lty=c(1,2))

}
\keyword{survey}

