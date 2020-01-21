library(survey)
library(AER)

load("cigsw.rda")

des<-svydesign(id=~1, weights=~wt, data=cigsw)
m<-svyivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi), design=des)

all.equal(as.vector(coef(m)), c(10.42009 , -1.588135, 0.6140887),tolerance=1e-6)
all.equal(as.vector(SE(m)), c( 1.047699, .3394232,  .3614382 ),tolerance=1e-6)
