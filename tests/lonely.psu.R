
## lonely PSUs by design
library(survey)
data(api)
## not certainty PSUs by fpc
ds<-svydesign(id = ~1, weights = ~pw, strata = ~dnum, data = apiclus1)
summary(ds)

options(survey.lonely.psu="fail")
try(svymean(~api00,ds))
try(svymean(~api00, as.svrepdesign(ds)))
options(survey.lonely.psu="remove")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="certainty")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="adjust")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="average")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))

## fpc specified
fpc<-ifelse(apiclus1$dnum==413, 1,1000)
ds<-svydesign(id = ~1, weights = ~pw, strata = ~dnum, data = apiclus1,fpc=fpc)
summary(ds)

options(survey.lonely.psu="fail")
try(svymean(~api00,ds))
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="remove")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="certainty")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="adjust")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))
options(survey.lonely.psu="average")
svymean(~api00,ds)
svymean(~api00, as.svrepdesign(ds))

rs<-as.svrepdesign(ds)
svytotal(~api00,rs)
SE(svytotal(~api00,subset(rs, dnum==413)))==0

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

## with adjustment
options(survey.adjust.domain.lonely=TRUE)
ds<-svydesign(id = ~1, weights = ~pw, strata = ~dnum, data = subset(apiclus1,dnum !=413))
ds1<-ds[-31,]
summary(ds1)

options(survey.lonely.psu="fail")
try(svymean(~api00,ds1))
options(survey.lonely.psu="remove")
svymean(~api00,ds1)
options(survey.lonely.psu="certainty")
svymean(~api00,ds1)
options(survey.lonely.psu="adjust")
svymean(~api00,ds1)
options(survey.lonely.psu="average")
svymean(~api00,ds1)

## checks for `svytotal()`

df_w_singleton <- data.frame(
  Stratum = c(1, 1, 2, 2, 3),
  PSU = c(1, 2, 3, 4, 5),
  Design_Weight = c(10.5, 10.5, 20, 20, 15),
  Sex = c("M", "F", "F", "M",  "M"),
  Age = c(30.5, 40.5, 35, 52, 44),
  Height = c(6.2, 5.0, 5.3, 5.7, 5.5)
)

design_w_singleton <- survey::svydesign(
  data = df_w_singleton,
  ids = ~ PSU, strata = ~ Stratum,
  weights = ~ Design_Weight
)

options("survey.lonely.psu" = "remove")

stopifnot(all.equal(
  target = 126625,
  current = as.numeric(
    vcov(svytotal(x = ~ Age, design = design_w_singleton))
  )
))

options("survey.lonely.psu" = "certainty")

stopifnot(all.equal(
  target = 126625,
  current = as.numeric(
    vcov(svytotal(x = ~ Age, design = design_w_singleton))
  )
))

options("survey.lonely.psu" = "adjust")

stopifnot(all.equal(
  target = 127579.8,
  current = as.numeric(
    vcov(svytotal(x = ~ Age, design = design_w_singleton))
  ),
  scale = 127579.8, tolerance = 0.000001
))
