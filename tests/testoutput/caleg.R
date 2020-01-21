##
## Calibration examples
##


## Example of calibration to first-stage clusters
library(survey)
data(api)

clusters<-table(apiclus2$dnum)
clusters<-clusters[clusters>1 & names(clusters)!="639"]
apiclus2a<-subset(apiclus2, dnum %in% as.numeric(names(clusters)))

dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2a)

popclusters<-subset(apipop, dnum %in% as.numeric(names(clusters)))

pop<-lapply(as.numeric(names(clusters)), function(cluster) {
  colSums(model.matrix(~api99, model.frame(~api99, subset(popclusters, dnum %in% cluster))))})

names(pop)<-names(clusters)

dclus2g<-calibrate(dclus2, ~api99, pop,stage=1)

svymean(~api99, dclus2)
svymean(~api99, dclus2g)

round(svyby(~api99, ~dnum, design=dclus2, svymean),4)

round(svyby(~api99, ~dnum, design=dclus2g, svymean),4)

## Averaging to first stage

dclus1<- svydesign(id = ~dnum, weights = ~pw, data = apiclus1, fpc = ~fpc)
pop<-colSums(cbind(1,apipop$enroll),na.rm=TRUE)

dclus1g<-calibrate(dclus1, ~enroll, pop, aggregate=1)

svytotal(~enroll,dclus1g)
svytotal(~api.stu,dclus1g)

#variation within clusters should be zero
all.equal(0, max(ave(weights(dclus1g),dclus1g$cluster,FUN=var),na.rm=TRUE))

##bounded weights
 dclus1g<-calibrate(dclus1, ~enroll, pop)
 range(weights(dclus1g)/weights(dclus1))
 dclus1gb<-calibrate(dclus1, ~enroll, pop, bounds=c(.6,1.5))
 range(weights(dclus1gb)/weights(dclus1))

## Ratio estimators
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
svytotal(~api.stu,dstrat)
common<-svyratio(~api.stu, ~enroll, dstrat, separate=FALSE)
total.enroll<-sum(apipop$enroll,na.rm=TRUE)
predict(common, total=total.enroll)
dstratg<-calibrate(dstrat,~enroll-1, total.enroll, variance=1)
svytotal(~api.stu, dstratg)

## postStratify vs calibrate in stratified sample (Ben French)
set.seed(17)
dat<-data.frame(y=rep(0:1,each=100),x=rnorm(200)+2*rep(0:1,each=100),
                z=rbinom(200,1,.2), fpc=rep(c(100,10000),each=100))
dat$w<-ifelse(dat$y,dat$z,1-dat$z)
popw<-data.frame(w=c("0","1"), Freq=c(2000,8000))
 des<-svydesign(id=~1,fpc=~fpc, data=dat,strata=~y)
postStratify(des,~w,popw)->dps
dcal<-calibrate(des,~factor(w), pop=c(10000,8000))

all.equal(SE(svymean(~x,dcal)),SE(svymean(~x,dps)))

## missing data in calibrated design
dps$variables$z[1]<-NA
summary(svyglm(y~z+x,design=dps,family=quasibinomial))

## Ratio estimator using the heteroskedasticity parameter (Daniel Oehm)
# should match the ratio estmate above
dstratgh <- calibrate(dstrat,~enroll-1, total.enroll, variance=apistrat$enroll)
svytotal(~api.stu, dstratgh)

## individual boundary constraints as multiplicative values (Daniel Oehm)
bnds <- list(
  lower = c(1, 1, rep(-Inf, nrow(apistrat)-2)), 
  upper = c(1, 1, rep(Inf, nrow(apistrat)-2))) # the first two weights will remain unchanged the others are free to move
lapply(bnds, head)
dstratg1<-calibrate(dstrat, ~enroll-1, total.enroll, bounds = bnds, variance=apistrat$enroll)
svytotal(~api.stu, dstratg1)
head(weights(dstrat))
head(weights(dstratg1))
all.equal(weights(dstrat)[1:2], weights(dstratg1)[1:2])

## individual boundary constraints as constant values (Daniel Oehm)
bnds <- list(
  lower = c(44.21, 44.21, rep(-Inf, nrow(apistrat)-2)), 
  upper = c(44.21, 44.21, rep(Inf, nrow(apistrat)-2))) # the first two weights will remain unchanged
lapply(bnds, head)
dstratg2<-calibrate(dstrat, ~enroll-1, total.enroll, bounds = bnds, bounds.const = TRUE, variance=apistrat$enroll)
svytotal(~api.stu, dstratg2)
head(weights(dstrat))
head(weights(dstratg2))
all.equal(round(weights(dstrat)[1:2], 8), round(weights(dstratg2)[1:2]), 8) # minor rounding error but all good

# sparse matrix support (Daniel Oehm)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
pop.totals<-c(`(Intercept)`=6194, stypeH=755, stypeM=1018)
dclus1g<-calibrate(dclus1, ~stype, pop.totals)
svymean(~api00, dclus1g)
svytotal(~enroll, dclus1g)

pop.totals<-c(`(Intercept)`=6194, stypeH=755, stypeM=1018)
dclus1g<-calibrate(dclus1, ~stype, pop.totals, sparse = TRUE)
svymean(~api00, dclus1g)
svytotal(~enroll, dclus1g)

pop.totals<-c(`(Intercept)`=6194, stypeH=755, stypeM=1018)
dclus1g<-calibrate(dclus1, ~stype, pop.totals, sparse = TRUE, calfun = "raking")
svymean(~api00, dclus1g)
svytotal(~enroll, dclus1g)
