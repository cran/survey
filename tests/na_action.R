## from Terry Therneau
library(survey)
load("naa.rda")

fit1e <- svyglm( pseudo ~ age34 + ccr5 + factor(times),  design= adata.s,na.action=na.exclude)
fit1o <- svyglm( pseudo ~ age34 + ccr5 + factor(times),  design= adata.s)
all.equal(coef(fit1e),coef(fit1o))
all.equal(vcov(fit1e),vcov(fit1o))
