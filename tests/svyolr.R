library(survey)

################################################################################
# Example from svyolr: runs OK
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
dclus1<-update(dclus1, mealcat=cut(meals,c(0,25,50,75,100)))

m<-svyolr(mealcat~avg.ed+mobility+stype, design=dclus1)
m

## Use regTermTest for testing multiple parameters
test<-regTermTest(m, ~avg.ed+stype, method="LRT")
################################################################################

################################################################################
# If we wrap everything into a function: error b/c it looks for design variable
# dclus1 in .GlobalEnv

foo <- function(x){
  dclus1<-svydesign(id=~dnum, weights=~pw, data=x, fpc=~fpc)
  dclus1<-update(dclus1, mealcat=cut(meals,c(0,25,50,75,100)))

  m<-svyolr(mealcat~avg.ed+mobility+stype, design=dclus1)
  ## Use regTermTest for testing multiple parameters
  regTermTest(m, ~avg.ed+stype, method="LRT")
}

# OK
foo(apiclus1)

# Clean-up everything but apiclus1 and foo
rm(list = setdiff(ls(),c('apiclus1','foo','test')))

# Error
test2<-foo(apiclus1)
################################################################################
all.equal(test,test2)

