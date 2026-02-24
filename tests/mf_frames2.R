##
## test multiframe results against results from the Frames2 package
##

library(survey)

data(phoneframes)
A_in_frames<-cbind(1, DatA$Domain=="ab")
B_in_frames<-cbind(DatB$Domain=="ba",1)

Bdes_pps<-svydesign(id=~1, fpc=~ProbB, data=DatB,pps=ppsmat(PiklB))
Ades_pps <-svydesign(id=~1, fpc=~ProbA,data=DatA,pps=ppsmat(PiklA))

## optimal constant (Hartley) weighting
mf_pps<-multiframe(list(Ades_pps,Bdes_pps),list(A_in_frames,B_in_frames),theta=0.7417399) 
t1<-svytotal(~Lei,mf_pps)
stopifnot(all.equal(as.vector(coef(t1)), 53259.86947))
stopifnot(all.equal(as.vector(vcov(t1)), 1652534, tol=1e-4))

## dividing by the expected number of selections (BKA or HH estimator)
Awts<-cbind(1/DatA$ProbA, ifelse(DatA$ProbB==0,0,1/DatA$ProbB))
Bwts<-cbind(ifelse(DatB$ProbA==0,0,1/DatB$ProbA),1/DatB$ProbB )

mf_pps2<-multiframe(list(Ades_pps,Bdes_pps),list(Awts,Bwts),estimator="expected") 
t2<-svytotal(~Lei,mf_pps2)
stopifnot(all.equal(as.vector(coef(t2)), 50953.07595))
stopifnot(all.equal(as.vector(vcov(t2)), 4116803, tol=1e-4))
