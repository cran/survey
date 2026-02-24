##
## internal multiframe checks of subpopulations
##


library(survey)

data(phoneframes)
A_in_frames<-cbind(1, DatA$Domain=="ab")
B_in_frames<-cbind(DatB$Domain=="ba",1)

Bdes_pps<-svydesign(id=~1, fpc=~ProbB, data=DatB,pps=ppsmat(PiklB))
Ades_pps <-svydesign(id=~1, fpc=~ProbA,data=DatA,pps=ppsmat(PiklA))
mf_pps<-multiframe(list(Ades_pps,Bdes_pps),list(A_in_frames,B_in_frames),theta=0.7417399) 

glm1<-svyglm(Lei~0+Domain, design=mf_pps)
m1a <- svymean(~Lei, subset(mf_pps, Domain=="a"))
stopifnot(all.equal(as.vector(coef(glm1))[1], as.vector(coef(m1a))))
stopifnot(all.equal(as.vector(SE(glm1)[1]), as.vector(SE(m1a))))


stopifnot(all.equal(coef(svymean(~Lei, subset(mf_pps, Domain=="a"))),
                    coef(svymean(~Lei, subset(Ades_pps, Domain=="a")))))

stopifnot(all.equal(as.vector(SE((svymean(~Lei, subset(mf_pps, Domain=="a"))))),
                    as.vector(SE(svymean(~Lei, subset(Ades_pps, Domain=="a"))))))


glm1<-svyglm(Lei~0+Domain, design=mf_pps)
m1a <- svymean(~Lei, subset(mf_pps, Domain=="a"))
all.equal(as.vector(coef(glm1))[1], as.vector(coef(m1a)))
all.equal(as.vector(SE(glm1)[1]), as.vector(SE(m1a)))


m2 <- svymean(~Lei, subset(mf_pps, Domain %in% c("a","ab")))
m2a<-svymean(~Lei, Ades_pps)

