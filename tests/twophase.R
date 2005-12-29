library(survey)

## two-phase simple random sampling.
data(pbc, package="survival")
pbc$id<-1:nrow(pbc)
(d2pbc<-twophase(id=list(~id,~id), data=pbc, subset=~I(trt==-9)))
m<-svymean(~bili, d2pbc)
all.equal(coef(m),with(pbc, mean(bili[trt==-9])))
all.equal(SE(m),
          with(pbc, sd(bili[trt==-9])/sqrt(sum(trt==-9))),
          tolerance=0.001)

## two-stage sampling as two-phase
data(mu284)
ii<-with(mu284, c(1:15, rep(1:5,n2[1:5]-3)))
mu284.1<-mu284[ii,]
mu284.1$id<-1:nrow(mu284.1)
mu284.1$sub<-rep(c(TRUE,FALSE),c(15,34-15))
dmu284<-svydesign(id=~id1+id2,fpc=~n1+n2, data=mu284)
## first phase cluster sample, second phase stratified within cluster
(d2mu284<-twophase(id=list(~id1,~id),strata=list(NULL,~id1),
                   fpc=list(~n1,NULL),data=mu284.1,subset=~sub))
summary(d2mu284)
t1<-svytotal(~y1, dmu284)
t2<-svytotal(~y1, d2mu284)
m1<-svymean(~y1, dmu284)
m2<-svymean(~y1, d2mu284)
all.equal(coef(t1),coef(t2))
all.equal(coef(m1),coef(m2))
all.equal(SE(t1),SE(t2))
all.equal(SE(t1),SE(t2))

## case-cohort design
##this example requires R 2.2.0 or later for cch and data.
library("survival")
data(nwtco, package="survival")
## unstratified, equivalent to Lin & Ying (1993)
(dcchs<-twophase(id=list(~seqno,~seqno), strata=list(NULL,~rel),
                 subset=~I(in.subcohort | rel), data=nwtco))
cch1<-svycoxph(Surv(edrel,rel)~factor(stage)+factor(histol)+I(age/12),
               design=dcchs)
## Using survival::cch 
subcoh <- nwtco$in.subcohort
selccoh <- with(nwtco, rel==1|subcoh==1)
ccoh.data <- nwtco[selccoh,]
ccoh.data$subcohort <- subcoh[selccoh]
cch2<-cch(Surv(edrel, rel) ~ factor(stage) + factor(histol) + I(age/12),
          data =ccoh.data, subcoh = ~subcohort, id=~seqno,
          cohort.size=4028, method="LinYing")

all.equal(coef(cch1),coef(cch2))
## variances equal only to just over 3 digits,
##  probably because of model-based variance in survival::cch
all.equal(as.vector(SE(cch1)),as.vector(SE(cch2)),tolerance=0.0005)

## bug report from Takahiro Tsuchiya for version 3.4
## We do not match Sarndal exactly, because our phase-one
## estimator has O(1/n.phase.2) bias.
rei<-read.table(textConnection(
"  id   N n.a h n.ah n.h   sub  y
1   1 300  20 1   12   5  TRUE  1
2   2 300  20 1   12   5  TRUE  2
3   3 300  20 1   12   5  TRUE  3
4   4 300  20 1   12   5  TRUE  4
5   5 300  20 1   12   5  TRUE  5
6   6 300  20 1   12   5 FALSE NA
7   7 300  20 1   12   5 FALSE NA
8   8 300  20 1   12   5 FALSE NA
9   9 300  20 1   12   5 FALSE NA
10 10 300  20 1   12   5 FALSE NA
11 11 300  20 1   12   5 FALSE NA
12 12 300  20 1   12   5 FALSE NA
13 13 300  20 2    8   3  TRUE  6
14 14 300  20 2    8   3  TRUE  7
15 15 300  20 2    8   3  TRUE  8
16 16 300  20 2    8   3 FALSE NA
17 17 300  20 2    8   3 FALSE NA
18 18 300  20 2    8   3 FALSE NA
19 19 300  20 2    8   3 FALSE NA
20 20 300  20 2    8   3 FALSE NA
"), header=TRUE)

des.rei <- twophase(id=list(~id,~id), strata=list(NULL,~h),
                    fpc=list(~N,NULL), subset=~sub, data=rei)
tot<- svytotal(~y, des.rei)

## based on Sarndal et al (9.4.14)
rei$w.ah <- rei$n.ah / rei$n.a
a.rei <- aggregate(rei, by=list(rei$h), mean, na.rm=TRUE)
a.rei$S.ysh <- tapply(rei$y, rei$h, var, na.rm=TRUE)
a.rei$y.u <- sum(a.rei$w.ah * a.rei$y)
V <- with(a.rei, sum(N * (N-1) * ((n.ah-1)/(n.a-1) - (n.h-1)/(N-1)) * w.ah * S.ysh / n.h))
V <- V + with(a.rei, sum(N * (N-n.a) * w.ah * (y - y.u)^2 / (n.a-1)))

a.rei$f.h<-with(a.rei, n.h/n.ah)
Vphase2<-with(a.rei, sum(N*N*w.ah^2* ((1-f.h)/n.h) *S.ysh))

a.rei$f<-with(a.rei, n.a/N)
a.rei$delta.h<-with(a.rei, (1/n.h)*(n.a-n.ah)/(n.a-1))
Vphase1<-with(a.rei, sum(N*N*((1-f)/n.a)*( w.ah*(1-delta.h)*S.ysh+ ((n.a)/(n.a-1))*w.ah*(y-y.u)^2)))

V
Vphase1
Vphase2
vcov(tot)
## phase 2 identical
all.equal(Vphase2,drop(attr(vcov(tot),"phases")$phase2))
## phase 1 differs by 2.6%
Vphase1/attr(vcov(tot),"phases")$phase1
