\name{api}
\alias{api}
\alias{apipop}
\alias{apiclus1}
\alias{apiclus2}
\alias{apistrat}
\docType{data}
\title{Student performance in California schools}
\description{
The Academic Performance Index is computed for all California schools
based on standardised testing of students. The data sets contain
information for all schools with at least 100 students and for various
probability samples of the data.
}
\usage{
data(apipop)
data(apistrat)
data(apiclus1)
data(apiclus2)
}
\format{
  The full population data in \code{pop} are a data frame with 6194 observations on the following 37 variables.
  \describe{
    \item{cds}{Unique identifier}
    \item{stype}{Elementary/Middle/High School}
    \item{name}{School name (15 characters)}
    \item{sname}{School name (40 characters)}
    \item{snum}{School number}
    \item{dname}{District name}
    \item{dnum}{District number}
    \item{cname}{County name}
    \item{cnum}{County number}
    \item{flag}{reason for missing data}
    \item{pcttest}{percentage of students tested}
    \item{api00}{API in 2000}
    \item{api99}{API in 1999}
    \item{target}{target for change in API}
    \item{growth}{Change in API}
    \item{sch.wide}{Met school-wide growth target?}
    \item{comp.imp}{Met Comparable Improvement target}
    \item{both}{Met both targets}
    \item{awards}{Eligible for awards program}
    \item{meals}{Percentage of students eligible for subsidized meals}
    \item{ell}{`English Language Learners' (percent)}
    \item{yr.rnd}{Year-round school}
    \item{mobility}{percentage of students for whom this is the first
      year at the school}
    \item{acs.k3}{average class size years K-3}
    \item{acs.46}{average class size years 4-6}
    \item{acs.core}{Number of core academic courses}
    \item{pct.resp}{percent where parental education level is known}
    \item{not.hsg}{percent parents not high-school graduates}
    \item{hsg}{percent parents who are high-school graduates}
    \item{some.col}{percent parents with some college}
    \item{col.grad}{percent parents with college degree}
    \item{grad.sch}{percent parents with postgraduate education}
    \item{avg.ed}{average parental education level}
    \item{full}{percent fully qualified teachers}
    \item{emer}{percent teachers with emergency qualifications}
    \item{enroll}{number of students enrolled}
    \item{api.stu}{number of students tested.}
  }
  The other data sets contain additional variables \code{pw} for
  sampling weights and \code{fpc} to compute finite population
  corrections to variance.
}
\details{
  \code{pop} is the entire population,
 \code{clus1} is a cluster sample of school districts, \code{strat} is
  a sample stratified by \code{stype}, and \code{clus2} is a two-stage
  cluster sample of schools within districts. 
}
\source{
  Data were obtained from the survey sampling help pages of UCLA
  Academic Technology Services, at \url{http://www.ats.ucla.edu/stat/stata/Library/svy_survey.htm}.
}
\references{
  The API program and original data files are at  \url{http://api.cde.ca.gov/}
}
\examples{
library(survey)
data(api)
mean(apipop$api00)
sum(apipop$enroll)

#stratified sample
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
summary(dstrat)
svymean(~api00, dstrat)
svytotal(~enroll, dstrat)

# one-stage cluster sample
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
summary(dclus1)
svymean(~api00, dclus1)
svytotal(~enroll, dclus1)

# two-stage cluster sample
dclus2<-svydesign(id=~dnum+snum, weights=~pw, data=apiclus2)
summary(dclus2)
svymean(~api00, dclus2)
svytotal(~enroll, dclus2)

# convert to replicate weights
rclus1<-as.svrepdesign(dclus1)
summary(rclus1)
svrepmean(~api00, rclus1)
svreptotal(~enroll, rclus1)

# post-stratify on school type
pop.types<-xtabs(~stype, data=apipop)

rclus1p<-postStratify(rclus1, ~stype, pop.types)
summary(rclus1p)
svrepmean(~api00, rclus1p)
svreptotal(~enroll, rclus1p)

}
\keyword{datasets}
