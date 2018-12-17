\name{svycralpha}
\alias{svycralpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cronbach's alpha
}
\description{
Compute Cronbach's alpha coefficient of reliability from survey data.  The formula is equation (2) of Cronbach (1951) only with design-based estimates of the variances.
}
\usage{
svycralpha(formula, design, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
One-sided formula giving the variables that make up the total score
}
  \item{design}{
survey design object
}
  \item{na.rm}{
\code{TRUE} to remove missing values
}
}

\value{
A number
}
\references{
Cronbach LJ (1951). "Coefficient alpha and the internal structure of tests". Psychometrika. 16 (3): 297-334. doi:10.1007/bf02310555.
}

\examples{
data(api)
dstrat<-svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat, 
    fpc = ~fpc)
svycralpha(~ell+mobility+avg.ed+emer+meals, dstrat)    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}% use one of  RShowDoc("KEYWORDS")

