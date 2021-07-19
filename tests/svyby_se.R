library(survey)
data(api)
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)

a<-svyby(~enroll,~stype,design=dstrat,svytotal,vartype=c("ci","se"))
b<-svyby(~enroll,~stype,design=dstrat,svytotal,vartype=c("se","ci"))


stopifnot(all.equal(SE(a),SE(b)))
