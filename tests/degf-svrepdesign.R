library(survey)
data(scd)

repweights<-2*cbind(c(1,0,1,0,1,0), c(1,0,0,1,0,1), c(0,1,1,0,0,1),
c(0,1,0,1,1,0))
scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights)

stopifnot(degf(scdrep)==3)

scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights, degf=4)
stopifnot(degf(scdrep)==4)

scdrep<-svrepdesign(data=scd, type="BRR", repweights=repweights, degf=2)
stopifnot(degf(scdrep)==2)

msg<-tryCatch(scdrep<-svrepdesign(data=scd, type="BRR",weights=~I(1000+0*ESA), repweights=repweights, combined.weights=FALSE,degf=10), 
   warning=function(w) w)

stopifnot(inherits(msg,"warning"))



