library(survey)
data(api)
dclus1 <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

dclus1<-update(dclus1, mealcat=cut(meals,c(0,25,50,75,100)))

## population marginal totals for each stratum
pop.types <- data.frame(stype=c("E","H","M"), Freq=c(4421,755,1018))
pop.schwide <- data.frame(sch.wide=c("No","Yes"), Freq=c(1072,5122))


## rake with the population totals
dclus1r<-rake(dclus1, list(~stype,~sch.wide), list(pop.types, pop.schwide))

# works
m <- svyolr(mealcat~avg.ed+mobility+stype, design=dclus1)

# fails in 4.1 (should work because svyolr's default na.action is na.omit)
m2 <- svyolr(mealcat~avg.ed+mobility+stype, design=dclus1r)

# fails in 4.1 (should work because NA values are subsetted out)
m3 <- svyolr(mealcat~avg.ed+mobility+stype, design=subset( dclus1r , !is.na( avg.ed ) ) )
