library(survey)
data(api)
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)
rclus2<-as.svrepdesign(dclus2)

m<-svyglm(I(comp.imp=="Yes")~1, design=dclus2, family=quasibinomial)
if(anyNA(confint(m, method="likelihood"))) stop("NA in confint")
mrep<-svyglm(I(comp.imp=="Yes")~1, design=rclus2, family=quasibinomial)
if(anyNA(confint(mrep))) stop("NA in confint")
if(anyNA(confint(mrep, method="likelihood"))) stop("NA in confint")

repstats <- withReplicates(rclus2, function(weights, data) {
  tot <- sum(data[['meals']] * weights)
  log_tot <- log(tot)
  return(c(tot, log_tot))
})
confint(repstats, level = 0.95)
