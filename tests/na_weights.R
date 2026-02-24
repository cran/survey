library(survey)

data(api)

dsrs <- svydesign(id=~1, weights=~pw, data=apisrs)

# Add records with 0 weight
apisrs_0 <- apisrs
apisrs_0$pw <- NA
apisrs2 <- rbind(apisrs, apisrs_0)
dsrs2 <- svydesign(id=~1, weights=~pw, data=apisrs2, na_weight="allow")

stopifnot(coef(svytotal(~enroll, dsrs))==coef(svytotal(~enroll, dsrs2)))
stopifnot(SE(svytotal(~enroll, dsrs))==SE(svytotal(~enroll, dsrs2)))
