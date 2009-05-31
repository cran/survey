library(survey)
data(election)

dpps<- svydesign(id=~1, weights=~wt, fpc=~p, data=election_pps, pps="brewer")
dppswr <-svydesign(id=~1, weights=~wt, data=election_pps)
svytotal(~Bush+Kerry+Nader, dpps)
svytotal(~Bush+Kerry+Nader, dppswr)

##subsets
svytotal(~Bush+Kerry+Nader, subset(dpps, Nader>0))

##multistage: should agree with STRS analysis
data(api)
dclus2<-svydesign(id = ~dnum + snum, fpc = ~fpc1 + fpc2, data = apiclus2)
dclus2pps<-svydesign(id = ~dnum + snum, fpc = ~I(40/fpc1) + I(pmin(1,5/fpc2)), data = apiclus2)

all.equal(svytotal(~sch.wide,dclus2), svytotal(~sch.wide,dclus2pps))
all.equal(svymean(~sch.wide,dclus2), svymean(~sch.wide,dclus2pps))
all.equal(svytotal(~enroll,dclus2), svytotal(~enroll,dclus2pps))

