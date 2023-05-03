
library(survey)
data(api)
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
a<-summary(svyglm(api00~ell+meals+mobility, design=dstrat))
b<- summary(svyivreg(api00~ell+meals+mobility, design=dstrat))
stopifnot(isTRUE(all.equal(a$cov.scaled, b$vcov)))
