\name{regTermTest}
\alias{regTermTest}
\alias{print.regTermTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wald test for a term in a regression model}
\description{
  Provides a Wald test of the hypothesis that all coefficients
  associated with a particular regression term are zero (or have some
  other specified values). Particularly useful as a substitute for
  \code{\link{anova}} when not fitting by maximum likelihood.
}
\usage{
regTermTest(model, test.terms, null=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A model object with \code{\link{coef}} and \code{\link{vcov}} methods}
  \item{test.terms}{Character string or one-sided formula giving name of
    term or terms to test}
  \item{null}{Null hypothesis values for parameters. Default is zeros}
}
\value{
  An object of class \code{regTermTest}.
}

\seealso{\code{\link{anova}}, \code{\link{vcov}}, \code{\link{contrasts}}}
\examples{
 data(esoph)
 model1 <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp * 
     alcgp, data = esoph, family = binomial())
 anova(model1)

 regTermTest(model1,"tobgp")
 regTermTest(model1,"tobgp:alcgp")
 regTermTest(model1, ~alcgp+tobgp:alcgp)
}
\keyword{regression}% at least one, from doc/KEYWORDS

