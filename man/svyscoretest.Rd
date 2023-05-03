\name{svyscoretest}
\alias{svyscoretest}
\alias{svyscoretest.svyglm}

\title{
Score tests in survey regression models
}
\description{
Performs two versions of the efficient score test. These are the same
for a single parameter. In the \code{working} score test, different
parameters are weighted according to the inverse of the estimated population Fisher
information. In the \code{pseudoscore} test, parameters are weighted according to the
inverse of their estimated covariance matrix.
}
\usage{
svyscoretest(model, drop.terms=NULL, add.terms=NULL,
method=c("working","pseudoscore","individual"),ddf=NULL,
lrt.approximation = "satterthwaite", ...)
\S3method{svyscoretest}{svyglm}(model, drop.terms=NULL, add.terms=NULL,
method=c("working","pseudoscore","individual"), ddf=NULL,
lrt.approximation = "satterthwaite",fullrank=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    A model of a class having a \code{svyscoretest} method (currently
    just \code{svyglm})
}
  \item{drop.terms}{
    Model formula giving terms to remove from \code{model}
  }
  \item{add.terms}{Model formula giving terms to add to \code{model} }
\item{method}{
  The type of score test to use. For a single parameter they are
  equivalent. To report tests for each column separately use \code{individual}
}
\item{ddf}{denominator degrees of freedom for an F or linear combination
  of F distributions. Use \code{Inf} to get chi-squared
  distributions. \code{NULL} asks for the model residual degrees of
  freedom, which is conservative. }
  \item{lrt.approximation}{
For the working score, the method for computing/approximating the null
distribution: see \code{\link{pchisqsum}}
  }
\item{fullrank}{If \code{FALSE} and \code{method="individual"}, keep
  even linearly dependent columns of the
  efficient score}
  \item{\dots}{
    for future expansion
}
}
\details{
The \code{working} score test will be asymptotically equivalent to the
Rao-Scott likelihood ratio test computed by \code{\link{regTermTest}}
and \code{\link{anova.svyglm}}. The paper by Rao, Scott and Skinner calls this
a "naive" score test. The null distribution is a linear combination
of chi-squared (or F) variables. 

The \code{pseudoscore} test will be
asymptotically equivalent to the Wald test computed by
\code{\link{regTermTest}}; it has a chi-squared (or F) null
distribution.

If \code{ddf} is negative or zero, which can happen with large numbers
of predictors and small numbers of PSUs, it will be changed to 1 with a warning.

}
\value{
For "pseudoscore" and "working" score methods, a named vector with the test
statistic, degrees of freedom, and p-value. For "individual" an object
of class "svystat"
}
\references{
JNK Rao, AJ Scott, and C
Rao, J., Scott, A., & Skinner, C. (1998). QUASI-SCORE TESTS WITH SURVEY
DATA. Statistica Sinica, 8(4), 1059-1070.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{regTermTest}}, \code{\link{anova.svyglm}}

}
\examples{
data(myco)
dmyco<-svydesign(id=~1, strata=~interaction(Age,leprosy),weights=~wt,data=myco)

m_full<-svyglm(leprosy~I((Age+7.5)^-2)+Scar, family=quasibinomial, design=dmyco)
svyscoretest(m_full, ~Scar)

svyscoretest(m_full,add.terms= ~I((Age+7.5)^-2):Scar)
svyscoretest(m_full,add.terms= ~factor(Age), method="pseudo")
svyscoretest(m_full,add.terms= ~factor(Age),method="individual",fullrank=FALSE)

svyscoretest(m_full,add.terms= ~factor(Age),method="individual")

}
