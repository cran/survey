
 library(survey)
example(svycoxph, ask=FALSE)
m<-update(model, .~.+I(protime^2))
a<-anova(m,model)
b<-anova(m, model,force=TRUE)
stopifnot(isTRUE(all.equal(b[2:6],a[c(3,4,6,7,8)])))
