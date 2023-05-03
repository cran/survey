\name{myco}
\alias{myco}
\docType{data}
\title{
Association between leprosy and BCG vaccination
}
\description{
These data are in a paper by JNK Rao and colleagues, on score tests
for complex survey data.  External information (not further specified) suggests
the functional form for the \code{Age} variable.
}
\usage{data("myco")}
\format{
  A data frame with 516 observations on the following 6 variables.
  \describe{
    \item{\code{Age}}{Age in years at the midpoint of six age strata}
    \item{\code{Scar}}{Presence of a BCG vaccination scar}
    \item{\code{n}}{Sampled number of cases (and thus controls) in the age stratum}
    \item{\code{Ncontrol}}{Number of non-cases in the population}
    \item{\code{wt}}{Sampling weight}
    \item{\code{leprosy}}{case status 0/1}
  }
}
\details{
The data are a simulated stratified case-control study drawn from a
population study conducted in a region of Malawi (Clayton and Hills,
1993, Table 18.1). The goal was to examine whether BCG vaccination against
tuberculosis protects against leprosy (the causative agents are both species of
_Mycobacterium_). Rao et al have a typographical error: the number of
non-cases in the population in the 25-30 age stratum is given as 4981
but 5981 matches both the computational output and the data as given by Clayton
and Hills.
}

\source{
JNK Rao, AJ Scott, and Skinner, C. (1998). QUASI-SCORE TESTS WITH SURVEY
DATA. Statistica Sinica, 8(4), 1059-1070.

Clayton, D., & Hills, M. (1993). Statistical Models in Epidemiology. OUP
}

\examples{
data(myco)
dmyco<-svydesign(id=~1, strata=~interaction(Age,leprosy),weights=~wt,data=myco)

m_full<-svyglm(leprosy~I((Age+7.5)^-2)+Scar, family=quasibinomial, design=dmyco)
m_age<-svyglm(leprosy~I((Age+7.5)^-2), family=quasibinomial, design=dmyco)
anova(m_full,m_age)

## unweighted model does not match
m_full
glm(leprosy~I((Age+7.5)^-2)+Scar, family=binomial, data=myco)

}
\keyword{datasets}
