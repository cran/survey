\name{poisson_sampling}
\alias{poisson_sampling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specify Poisson sampling design
}
\description{
Specify a design where units are sampled independently from the population, with known probabilities. This design is often used theoretically, but is rarely used in practice because the sample size is variable. This function calls \code{\link{ppscov}} to specify a sparse sampling covariance matrix. 
}
\usage{
poisson_sampling(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
Vector of sampling probabilities
}
}

\value{
Object of class \code{ppsdcheck}
}

\seealso{
\code{\link{ppscov}}, \code{\link{svydesign}}
}
\examples{
data(api)
apipop$prob<-with(apipop, 200*api00/sum(api00))
insample<-rbinom(nrow(apipop),1,apipop$prob)
apipois<-apipop[insample,]
despois<-svydesign(id=~1, prob=~prob, pps=poisson_sampling(apipois$prob), data=apipois)

svytotal(~api00, despois)

## SE formula
sqrt(sum( (apipois$api00*weights(despois))^2*(1-apipois$prob)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{survey }
