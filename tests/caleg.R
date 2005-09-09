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
