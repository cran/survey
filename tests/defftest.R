library(survey)

data(api)

## one-stage cluster sample
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

# svyglm model
mod <- svyglm(api99 ~ enroll + api.stu, design = dclus1, deff = TRUE)

#deffs returned from svyglm model - implausibly high
deff(mod)
#> (Intercept)      enroll     api.stu 
#>    351.3500    457.6799    491.5567

# run mod with same data and glm()
srs_mod <- glm(api99 ~ enroll + api.stu, data = apiclus1)

# manually calculate deffs

clust_se <- summary(mod)$coefficients[,2]
srs_se <- summary(srs_mod)$coefficients[,2]

deffs <- clust_se^2 / srs_se^2
stopifnot(all.equal(deffs, deff(mod)))

