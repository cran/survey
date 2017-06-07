\name{api}
\alias{api}
\alias{apipop}
\alias{apiclus1}
\alias{apiclus2}
\alias{apistrat}
\alias{apisrs}
\docType{data}
\title{Student performance in California schools}
\description{
The Academic Performance Index is computed for all California schools
based on standardised testing of students. The data sets contain
information for all schools with at least 100 students and for various
probability samples of the data.
}
\usage{
data(api)
}
\format{
  The full population data in \code{apipop} are a data frame with 6194 observations on the following 37 variables.
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
  \code{apipop} is the entire population, \code{apisrs} is a simple random sample,
 \code{apiclus1} is a cluster sample of school districts, \code{apistrat} is
  a sample stratified by \code{stype}, and \code{apiclus2} is a two-stage
  cluster sample of schools within districts.  The sampling weights in
 \code{apiclus1} are incorrect (the weight should be 757/15) but are as 
  obtained from UCLA.
}
\source{
  Data were obtained from the survey sampling help pages of UCLA
  Academic Technology Services; these pages are no longer on line. 
}
\references{
  The API program and original data files are at  \url{http://api.cde.ca.gov/}
}
\examples{
library(survey)
data(api)
mean(apipop$api00)
sum(apipop$enroll, na.rm=TRUE)

#stratified sample
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
summary(dstrat)
svymean(~api00, dstrat)
svytotal(~enroll, dstrat, na.rm=TRUE)

# one-stage cluster sample
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
summary(dclus1)
svymean(~api00, dclus1)
svytotal(~enroll, dclus1, na.rm=TRUE)

# two-stage cluster sample
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)
summary(dclus2)
svymean(~api00, dclus2)
svytotal(~enroll, dclus2, na.rm=TRUE)

# two-stage `with replacement'
dclus2wr<-svydesign(id=~dnum+snum, weights=~pw, data=apiclus2)
summary(dclus2wr)
svymean(~api00, dclus2wr)
svytotal(~enroll, dclus2wr, na.rm=TRUE)


# convert to replicate weights
rclus1<-as.svrepdesign(dclus1)
summary(rclus1)
svymean(~api00, rclus1)
svytotal(~enroll, rclus1, na.rm=TRUE)

# post-stratify on school type
pop.types<-xtabs(~stype, data=apipop)

rclus1p<-postStratify(rclus1, ~stype, pop.types)
dclus1p<-postStratify(dclus1, ~stype, pop.types)
summary(dclus1p)
summary(rclus1p)

svymean(~api00, dclus1p)
svytotal(~enroll, dclus1p, na.rm=TRUE)

svymean(~api00, rclus1p)
svytotal(~enroll, rclus1p, na.rm=TRUE)

}
\keyword{datasets}
