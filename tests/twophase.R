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
