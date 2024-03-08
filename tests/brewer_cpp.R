## pps="brewer" can't use rcpp
## this checks that it doesn't
library(survey)
data(election)

dpps_br<- svydesign(id=~1,  fpc=~p, data=election_pps, pps="brewer")
options(survey.use_rcpp=TRUE)
a<-svytotal(~Bush+Kerry+Nader, dpps_br)
options(survey.use_rcpp=FALSE)
b<-svytotal(~Bush+Kerry+Nader, dpps_br)

stopifnot(identical(a,b))
