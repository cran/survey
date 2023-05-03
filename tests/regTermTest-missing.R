library(survey)
data(nhanes)
design <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE,data=nhanes)

a<-svyglm(formula = I(race == 1) ~ HI_CHOL + agecat + RIAGENDR,  design = subset(design,!is.na(HI_CHOL)), family=quasibinomial)
b<-svyglm(formula = I(race == 1) ~ HI_CHOL + agecat + RIAGENDR,  design =design , family=quasibinomial)

ta<-regTermTest(a, ~HI_CHOL)
tb<-regTermTest(b, ~HI_CHOL)

stopifnot(isTRUE(all.equal(ta$chisq, tb$chisq)))
stopifnot(isTRUE(all.equal(ta$lambda, tb$lambda)))
