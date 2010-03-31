\name{stratsample}
\Rdversion{1.1}
\alias{stratsample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Take a stratified sample
}
\description{
This function takes a stratified sample without replacement from a data set. 
}
\usage{
stratsample(strata, counts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strata}{
Vector of stratum identifiers; will be coerced to character
}
  \item{counts}{
  named vector of stratum sample sizes, with names corresponding to the values of \code{as.character(strata)}
}
}

\value{
vector of indices into \code{strata} giving the sample
}


\seealso{
\code{\link{sample}}

The "sampling" package has many more sampling algorithms.
}
\examples{
  data(api)
  s<-stratsample(apipop$stype, c("E"=5,"H"=4,"M"=2))
  table(apipop$stype[s])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

