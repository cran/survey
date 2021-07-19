## test use of replicates in svyby, svycontrast
library(survey)

data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1<-as.svrepdesign(dclus1)

meanlogs_without<-svyby(~log(enroll),~stype,svymean, design=rclus1,covmat=TRUE)
c_without<-svycontrast(meanlogs_without, quote(exp(E-H)))
vcov(c_without)

meanlogs_with<-svyby(~log(enroll),~stype,svymean, design=rclus1,covmat=TRUE,return.replicates=TRUE)

c_with<-svycontrast(meanlogs_with, quote(exp(E-H)))

v_with<- vcov(rclus1, c_with$replicates)

r<- attr(meanlogs_with, "replicates")
vr_with<-vcov(rclus1,exp(r[,1]-r[,2]))

stopifnot(all.equal(as.numeric(v_with),as.numeric(vr_with)))
stopifnot(all.equal(as.numeric(v_with),as.numeric(vcov(c_with))))
