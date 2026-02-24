
library(survey)
data(nwtco)
dcchs<-twophase(id = list(~seqno, ~seqno), strata = list(NULL, ~rel), 
    subset = ~I(in.subcohort | rel), data = nwtco)
mcchs<-multiphase(id = list(~seqno, ~seqno), strata = list(NULL, ~rel), 
    subset = list(~I(in.subcohort | rel)), probs = list(~1, NULL), 
    data = nwtco)
dcchs
mcchs

old<-svymean(~edrel, dcchs)
new<-svymean(~edrel, mcchs)
stopifnot(all.equal(coef(old), coef(new)))
stopifnot(all.equal(as.vector(SE(old)), as.vector(SE(new))))
stopifnot(all.equal( attr(vcov(old),"phases")[[1]], attr(vcov(new),"phases")[[1]]))
stopifnot(all.equal( attr(vcov(old),"phases")[[2]], attr(vcov(new),"phases")[[2]]))

old<-svytotal(~edrel, dcchs)
new<-svytotal(~edrel, mcchs)
stopifnot(all.equal(coef(old), coef(new)))
stopifnot(all.equal(as.vector(SE(old)), as.vector(SE(new))))
stopifnot(all.equal( attr(vcov(old),"phases")[[1]], attr(vcov(new),"phases")[[1]]))
stopifnot(all.equal( attr(vcov(old),"phases")[[2]], attr(vcov(new),"phases")[[2]]))

old<-svytotal(~rel, dcchs) ##stratification variable, has no phase-two variance
new<-svytotal(~rel, mcchs)
stopifnot(all.equal(coef(old), coef(new)))
stopifnot(all.equal(as.vector(SE(old)), as.vector(SE(new))))
stopifnot(all.equal( attr(vcov(old),"phases")[[1]], attr(vcov(new),"phases")[[1]]))
stopifnot(all.equal( attr(vcov(old),"phases")[[2]], attr(vcov(new),"phases")[[2]]))

## not identical yet, need to update
m<-calibrate(mcchs,~factor(stage)+rel, phase=2, calfun="raking")
d<-calibrate(dcchs,~factor(stage)+rel, phase=2, calfun="raking")
old<-vcov(svytotal(~factor(stage), m))
new<-vcov(svytotal(~factor(stage), d))
stopifnot(all(abs(attr(old,"phases")[[2]])<1e-10))
stopifnot(all(abs(attr(new,"phases")[[2]])<1e-10))
