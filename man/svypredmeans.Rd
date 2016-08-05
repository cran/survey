\name{svypredmeans}
\alias{svypredmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictive marginal means
}
\description{
Predictive marginal means for a generalised linear model, using the method of Korn and Graubard (1999) and matching the results of SUDAAN. The predictive marginal mean for one level of a factor is the probability-weighted average of the fitted values for the model on new data where all the observations are set to that level of the factor but have whatever values of adjustment variables they really have.
}
\usage{
svypredmeans(adjustmodel, groupfactor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjustmodel}{
A generalised linear model fit by \code{\link{svyglm}} with the adjustment variable but without the factor for which predictive means are wanted
}
  \item{groupfactor}{
A one-sided formula specifying the factor for which predictive means are wanted. Can use, eg, \code{~interaction(race,sex)} for combining variables}

}

\value{
An object of class \code{svystat} with the predictive marginal means and their covariance matrix. 
}
\references{
Graubard B, Korn E (1999) "Predictive Margins with Survey Data" Biometrics 55:652-659

Bieler, Brown, Williams, & Brogan (2010) "Estimating Model-Adjusted Risks, Risk Differences, and Risk Ratios From Complex Survey Data" Am J Epi DOI: 10.1093/aje/kwp440}
\note{
It is possible to supply an adjustment model with only an intercept, but the results are then the same as \code{\link{svymean}}

It makes no sense to have a variable in the adjustment model that is part of the grouping factor, and will give an error message or \code{NA}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{svyglm}}

Worked example using National Health Interview Survey data: \url{https://gist.github.com/tslumley/2e74cd0ac12a671d2724}
}
\examples{
data(nhanes)
nhanes_design <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE,data=nhanes)
agesexmodel<-svyglm(HI_CHOL~agecat+RIAGENDR, design=nhanes_design,family=quasibinomial)
## high cholesterol by race/ethnicity, adjusted for demographic differences
means<-svypredmeans(agesexmodel, ~race)
means
## relative risks compared to non-Hispanic white
svycontrast(means,quote(`1`/`2`))
svycontrast(means,quote(`3`/`2`))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{regression}% __ONLY ONE__ keyword per line
