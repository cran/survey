
library(survey)
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1<-as.svrepdesign(dclus1)

a<-svyby(~api00+api99, ~comp.imp+sch.wide,design=rclus1,svymean,
         covmat=TRUE,drop.empty.groups=FALSE)
b<-svyby(~api00+api99, ~comp.imp+sch.wide,design=rclus1,svymean,
         covmat=TRUE,drop.empty.groups=TRUE)

stopifnot(all(na.omit(
              as.vector(as.matrix(SE(a)))==sqrt(diag(vcov(a)))
)))
stopifnot(all(
              as.vector(as.matrix(SE(b)))==sqrt(diag(vcov(b)))
              ))
