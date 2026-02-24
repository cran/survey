## bug report from Thomas Leeper, fixed in version 3.32-3

library("survey")
data(api)
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)

# pass `family` directly (WORKS!)
svyglm(api00~ell+meals+mobility, design=dstrat, family = gaussian())

# passing `family` via ... (WORKS!)
myfun1 <- function(formula, design, ...) {
    svyglm(formula, design = design, ...)
}
myfun1(api00~ell+meals+mobility, design=dstrat, family = gaussian())

# passing `family` via default argument (DOES NOT WORK!)
myfun2 <- function(formula, design, family = gaussian()) {
    svyglm(formula, design = design, family = family)
}
myfun2(api00~ell+meals+mobility, design=dstrat, family = gaussian())
