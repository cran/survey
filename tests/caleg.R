##
## Example of calibration to first-stage clusters
##

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
