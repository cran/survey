\name{paley}
\alias{paley}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Paley-type Hadamard matrices}
\description{
  Computes a Hadamard matrix of dimension p+1, where p is a prime,
  and p+1 is a multiple of 4, using the Paley construction. Works
  for \code{n} up to 1367.
}
\usage{
paley(n, nmax = 2 * n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Minimum size for matrix}
  \item{nmax}{Maximum size for matrix}
}

\value{
  A matrix of zeros and ones, or \code{NULL} if no matrix smaller than
  \code{nmax} can be found.
}

\seealso{ \code{\link{hadamard}}}
\examples{

p<-paley(11)

p1<-2*p-1
## HH^T is diagonal for any Hadamard matrix
p1\%*\%t(p1)

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{algebra}

