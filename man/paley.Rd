\name{paley}
\alias{paley}
\alias{is.hadamard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Paley-type Hadamard matrices}
\description{
  Computes a Hadamard matrix of dimension \eqn{(p+1)\times 2^k}{(p+1)*2^k}, where p is a prime,
  and p+1 is a multiple of 4, using the Paley construction. Used by \code{\link{hadamard}}.
}
\usage{
paley(n, nmax = 2 * n, prime=NULL, check=!is.null(prime))

is.hadamard(H, style=c("0/1","+-"), full.orthogonal.balance=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Minimum size for matrix}
  \item{nmax}{Maximum size for matrix. Ignored if \code{prime} is specified.}
  \item{prime}{Optional. A prime at least as large as
    \code{n},  such that \code{prime+1} is divisible by 4.}
  \item{check}{Check that the resulting matrix is of Hadamard type}
  \item{H}{Matrix}
  \item{style}{\code{"0/1"} for a matrix of 0s and 1s, \code{"+-"} for a
    matrix of \eqn{\pm 1}{+/-1}.}
  \item{full.orthogonal.balance}{Require full orthogonal balance?}
}

\value{
  For \code{paley}, a matrix of zeros and ones, or \code{NULL} if no matrix smaller than
  \code{nmax} can be found.

  For \code{is.hadamard}, \code{TRUE} if \code{H} is a Hadamard matrix.
}
\details{
  The Paley construction gives a Hadamard matrix of order p+1 if p is
  prime and p+1 is a multiple of 4.  This is then expanded to order
  \eqn{(p+1)\times 2^k}{(p+1)*2^k} using the Sylvester construction.
  
  \code{paley} knows primes up to 7919.  The user can specify a prime
  with the \code{prime} argument, in which case a matrix of order
  \eqn{p+1}{p+1} is constructed.

  If \code{check=TRUE} the code uses \code{is.hadamard} to check that
  the resulting matrix really is of Hadamard type, in the same way as in
  the example below. As this test takes \eqn{n^3}{n^3} time it is
  preferable to just be sure that \code{prime} really is prime.

  A Hadamard matrix including a row of 1s gives BRR designs where the
  average of the replicates for a linear statistic is exactly the full
  sample estimate. This property is called full orthogonal balance.
}
\references{
Cameron PJ (2005) Hadamard Matrices. In: The
Encyclopedia of Design Theory \url{http://www.maths.qmul.ac.uk/~lsoicher/designtheory.org/library/encyc/}

}
\seealso{ \code{\link{hadamard}}}
\examples{

M<-paley(11)

is.hadamard(M)
## internals of is.hadamard(M)
H<-2*M-1
## HH^T is diagonal for any Hadamard matrix
H\%*\%t(H)

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{algebra}

