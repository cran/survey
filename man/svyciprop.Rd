\name{svyciprop}
\alias{svyciprop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence intervals for proportions }
\description{
Computes confidence intervals for proportions using methods that may be
more accurate near 0 and 1 than simply using \code{confint(svymean())}.
}
\usage{
svyciprop(formula, design, method = c("logit", "likelihood", "asin", "beta",
"mean","xlogit"), level = 0.95, df=degf(design),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula specifying a single binary variable}
  \item{design}{ survey design object}
  \item{method}{ See Details below. Partial matching is done on the argument.}
  \item{level}{Confidence level for interval}
  \item{df}{denominator degrees of freedom, for all methods except
    \code{"beta"}.  Use \code{Inf} for confidence intervals based on a
    Normal distribution, and for \code{"likelihood"} and \code{"logit"}
    use \code{NULL} for the default method in glms (currently
    \code{degf(design)-1}, but this may be improved in the future)}
  \item{\dots}{For \code{"mean"} and \code{"asin"}, this is passed to \code{\link{confint.svystat}} }
}
\details{
The \code{"logit"} method fits a logistic regression model and computes a
Wald-type interval on the log-odds scale, which is then transformed to
the probability scale.

The \code{"likelihood"} method uses the (Rao-Scott) scaled chi-squared distribution
for the loglikelihood from a binomial distribution.

The \code{"asin"} method uses the variance-stabilising transformation
for the binomial distribution, the arcsine square root, and then
back-transforms the interval to the probability scale

The \code{"beta"} method uses the incomplete beta function as in
\code{\link{binom.test}}, with an effective sample size based on the
estimated variance of the proportion. (Korn and Graubard, 1998)

The \code{"xlogit"} method uses a logit transformation of the mean and
then back-transforms to the probablity scale. This appears to be the
method used by SUDAAN and SPSS COMPLEX SAMPLES. 

The \code{"mean"} method is a Wald-type interval on the probability
scale, the same as \code{confint(svymean())}

All methods undercover for probabilities close enough to zero or one,
but \code{"beta"}, \code{"likelihood"}, \code{"logit"}, and \code{"logit"} are noticeably
better than the other two. None of the methods will work when the
observed proportion is exactly 0 or 1.

The \code{confint} method extracts the confidence interval; the
\code{vcov} and \code{SE} methods just report the variance or standard
error of the mean. 
}
\value{
 The point estimate of the proportion, with the confidence interval as
 an attribute
}
\references{
Rao, JNK, Scott, AJ (1984) "On Chi-squared Tests For Multiway
Contingency Tables with Proportions Estimated From Survey Data" Annals
of Statistics 12:46-60.

Korn EL, Graubard BI. (1998) Confidence Intervals For Proportions With
  Small Expected Number of Positive Counts Estimated From Survey
  Data. Survey Methodology 23:193-201.
  
}

\seealso{ \code{\link{svymean}}, \code{\link{yrbs}} }
\examples{
data(api)
dclus1<-svydesign(id=~dnum, fpc=~fpc, data=apiclus1)

svyciprop(~I(ell==0), dclus1, method="li")
svyciprop(~I(ell==0), dclus1, method="lo")
svyciprop(~I(ell==0), dclus1, method="as")
svyciprop(~I(ell==0), dclus1, method="be")
svyciprop(~I(ell==0), dclus1, method="me")
svyciprop(~I(ell==0), dclus1, method="xl")

## reproduces Stata svy: mean
svyciprop(~I(ell==0), dclus1, method="me", df=degf(dclus1))
## reproduces Stata svy: prop
svyciprop(~I(ell==0), dclus1, method="lo", df=degf(dclus1))


rclus1<-as.svrepdesign(dclus1)
svyciprop(~I(emer==0), rclus1, method="li")
svyciprop(~I(emer==0), rclus1, method="lo")
svyciprop(~I(emer==0), rclus1, method="as")
svyciprop(~I(emer==0), rclus1, method="be")
svyciprop(~I(emer==0), rclus1, method="me")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

