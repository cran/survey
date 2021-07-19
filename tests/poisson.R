## check poisson sampling
library(survey)
data(api)
set.seed(2021-7-15)
apipop$prob<-apipop$api00/1000
insample<-rbinom(nrow(apipop),1,apipop$prob)
apipois<-apipop[insample,]
des<-svydesign(id=~1, prob=~prob, pps=poisson_sampling(apipois$prob), data=apipois)

stopifnot(isTRUE(all.equal(
 as.vector(SE(svytotal(~api00,design=des))),
 as.vector(sqrt(sum( (apipois$api00*weights(des))^2*(1-apipois$prob))))
 )))
 