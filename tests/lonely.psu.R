
## lonely PSUs by design
library(survey)
data(api)
## not certainty PSUs by fpc
ds<-svydesign(id = ~1, weights = ~pw, strata = ~dnum, data = apiclus1)
summary(ds)

options(survey.lonely.psu="fail")
try(svymean(~api00,ds))
options(survey.lonely.psu="remove")
svymean(~api00,ds)
options(survey.lonely.psu="certainty")
svymean(~api00,ds)
options(survey.lonely.psu="adjust")
svymean(~api00,ds)
options(survey.lonely.psu="average")
svymean(~api00,ds)

## fpc specified
fpc<-ifelse(apiclus1$dnum==413, 1,1000)
ds<-svydesign(id = ~1, weights = ~pw, strata = ~dnum, data = apiclus1,fpc=fpc)
summary(ds)

options(survey.lonely.psu="fail")
try(svymean(~api00,ds))
options(survey.lonely.psu="remove")
svymean(~api00,ds)
options(survey.lonely.psu="certainty")
svymean(~api00,ds)
options(survey.lonely.psu="adjust")
svymean(~api00,ds)
options(survey.lonely.psu="average")
svymean(~api00,ds)

## lonely PSUs after subsetting
ds<-svydesign(id = ~1, weights = ~pw, strata = ~dnum, data = subset(apiclus1,dnum !=413))
ds1<-ds[-31,]
summary(ds1)

options(survey.lonely.psu="fail")
svymean(~api00,ds1)
options(survey.lonely.psu="remove")
svymean(~api00,ds1)
options(survey.lonely.psu="certainty")
svymean(~api00,ds1)
options(survey.lonely.psu="adjust")
svymean(~api00,ds1)
options(survey.lonely.psu="average")
svymean(~api00,ds1)
