##
## Domain means can be written as ratio estimators or as regression coefficients
##
## This code checks that subsetting the design object gives the same results as
## these approaches.
##


library(survey)
data(fpc)
dfpc<-svydesign(id=~psuid,strat=~stratid,weight=~weight,data=fpc,nest=TRUE)
dsub<-subset(dfpc,x>4)
(m1<-svymean(~x,design=dsub))

## These should give the same domain estimates and standard errors
(m2<-svyby(~x,~I(x>4),design=dfpc, svymean,keep.var=TRUE))
m3<-svyglm(x~I(x>4)+0,design=dfpc)
summary(m3)
(m4<-svyratio(~I(x*(x>4)),~as.numeric(x>4), dfpc))
stopifnot(isTRUE(all.equal(SE(m2), as.vector(SE(m3)))))
stopifnot(isTRUE(all.equal(SE(m2)[2], as.vector(SE(m4)))))

## with strata
data(api)
dstrat<-svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
m1<-svymean(~enroll, subset(dstrat, comp.imp=="Yes"))
m2<-svyglm(enroll~comp.imp-1, dstrat)
m3<- svyratio(~I(enroll*(comp.imp=="Yes")), ~as.numeric(comp.imp=="Yes"), dstrat)
stopifnot(isTRUE(all.equal(as.vector(SE(m2)["comp.impYes"]), as.vector(SE(m1)))))
stopifnot(isTRUE( all.equal(as.vector(SE(m1)), as.vector(drop(SE(m3))))))

## with calibration
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
pop.totals<-c(`(Intercept)`=6194, stypeH=755, stypeM=1018)
(dclus1g3 <- calibrate(dclus1, ~stype+api99, c(pop.totals, api99=3914069)))

m1<-svymean(~api00, subset(dclus1g3, comp.imp=="Yes"))
m3<-svyratio(~I(api00*(comp.imp=="Yes")), ~as.numeric(comp.imp=="Yes"), dclus1g3)
m2<-svyglm(api00~comp.imp-1, dclus1g3)
stopifnot(isTRUE( all.equal(as.vector(SE(m2)["comp.impYes"]), as.vector(SE(m1)))))
stopifnot(isTRUE( all.equal(as.vector(SE(m1)), as.vector(drop(SE(m3))))))

## with raking
pop.types <- data.frame(stype=c("E","H","M"), Freq=c(4421,755,1018))
pop.schwide <- data.frame(sch.wide=c("No","Yes"), Freq=c(1072,5122))
dclus1r<-rake(dclus1, list(~stype,~sch.wide), list(pop.types, pop.schwide))
m1<-svymean(~api00, subset(dclus1r, comp.imp=="Yes"))
m2<-svyglm(api00~comp.imp-1, dclus1r)
m3<-svyratio(~I(api00*(comp.imp=="Yes")), ~as.numeric(comp.imp=="Yes"), dclus1r)
stopifnot(isTRUE( all.equal(as.vector(SE(m2)["comp.impYes"]), as.vector(SE(m1)))))
stopifnot(isTRUE( all.equal(as.vector(SE(m1)), as.vector(drop(SE(m3))))))
