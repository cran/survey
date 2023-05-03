### R code from vignette source 'precalibrated.Rnw'

###################################################
### code chunk number 1: precalibrated.Rnw:16-19
###################################################
library(survey)
data(api)
dclus1 <- svydesign(id = ~dnum, weights = ~pw, data = apiclus1, fpc = ~fpc)


###################################################
### code chunk number 2: precalibrated.Rnw:24-28
###################################################
sum(weights(dclus1))
dim(apipop)
dclus1<-update(dclus1, one=rep(1,nrow(dclus1)))
svytotal(~one,dclus1)


###################################################
### code chunk number 3: precalibrated.Rnw:34-36
###################################################
cal_dclus1<-calibrate(dclus1, formula=~1, population=sum(weights(dclus1)))
svytotal(~one,cal_dclus1)


###################################################
### code chunk number 4: precalibrated.Rnw:40-41
###################################################
summary(weights(cal_dclus1)/weights(dclus1))


###################################################
### code chunk number 5: precalibrated.Rnw:45-50
###################################################
precal_dclus1<-svydesign(id = ~dnum, weights = ~pw, data = apiclus1,
                         fpc = ~fpc, calibrate.formula=~1)
precal_dclus1<-update(precal_dclus1, one=rep(1,nrow(dclus1)))

svytotal(~one,precal_dclus1)


###################################################
### code chunk number 6: precalibrated.Rnw:55-64
###################################################
(enroll_t<-svytotal(~enroll, dclus1))
(enroll_m<-svymean(~enroll, dclus1))
SE(enroll_m)
SE(enroll_t)/6194

(cenroll_t<-svytotal(~enroll, precal_dclus1))
(cenroll_m<-svymean(~enroll, precal_dclus1))
SE(cenroll_m)
SE(cenroll_t)/6194


