library(survey)
data(fpc)
dfpc<-svydesign(id=~psuid,strat=~stratid,weight=~weight,data=fpc,nest=TRUE)
dsub<-subset(dfpc,x>4)
(m1<-svymean(~x,design=dsub))

## These should give the same domain estimates and standard errors
(m2<-svyby(~x,~I(x>4),design=dfpc, svymean,keep.var=TRUE))
m3<-svyglm(x~I(x>4)+0,design=dfpc)
summary(m3)

stopifnot(identical(TRUE,all.equal(m2$SE, SE(m3))))

