\name{svyivreg}
\alias{svyivreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two-stage least-squares for instrumental variable regression
}
\description{
Estimates regressions with endogenous covariates using two-stage least squares. The function uses \code{ivreg} from the \code{AER} package for the main computations, and follows the syntax of that function.

}
\usage{
svyivreg(formula, design, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{formula specification(s) of the regression
          relationship and the instruments. See Details for details}
  \item{design}{
A survey design object
}
  \item{\dots}{
For future expansion
}
}
\details{
Regressors and instruments for \code{svyivreg} are specified
     in a formula with two parts on the right-hand side, e.g., \code{y ~ x1
     + x2 | z1 + z2 + z3}, where \code{x1} and \code{x2} are the regressors and
     \code{z1}, \code{z2}, and \code{z3} are the instruments. Note that exogenous
     regressors have to be included as instruments for themselves. For
     example, if there is one exogenous regressor \code{ex} and one
     endogenous regressor \code{en} with instrument \code{in}, the appropriate
     formula would be \code{y ~ ex + en | ex + in}. Equivalently, this can
     be specified as \code{y ~ ex + en | . - en + in}, i.e., by providing an
     update formula with a \code{.} in the second part of the formula. }
\value{
An object of class \code{svyivreg}
}
\references{
\url{https://notstatschat.rbind.io/2019/07/16/adding-new-functions-to-the-survey-package/}
}


\seealso{
\code{\link[AER]{ivreg}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey }% use one of  RShowDoc("KEYWORDS")
