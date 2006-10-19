##
##  standard errors using pairwise sampling probabilities
##

##
## don't need fpc, do need probabilities
##
## prob is either a list of lists of matrices, a list of Matrices,
##      a single matrix or a formula
##
##      if a formula, pairwise= specifies the approximation to use
##
##      if a single matrix or list of Matrices, don't want strata
##      if a single matrix, don't even want id.
##
## we compute and store \check{Delta}= (pi_ij/pi_ipi_j -1)/pi_ij
##  (and perhaps pi_ij)?

pairdesign<-function(id, prob, strata, data, pairwise=c("srs","overton")){
    
    if (is.formula(id))
        id<-model.frame(id, data)


        


}

##
## here prob is a vector, dcheck is a matrix, list of Matrices,
##                               or list of lists of matrices
##

pi2Dcheck<-function(pmat){
    (pmat-outer(diag(pmat)))/pmat
}

strat2Dcheck<-function(strata, prob){
    n<-length(strata)
    rval<-matrix(0, n,n)
    sampsize<-ave(strata,strata,FUN=length)
    strats<-unique(strata)
    for(strat in strats){
        these <- strata == strat
        rval[these,these]<- -(1-prob[these])/(sampsize[these]-1)
    }
    diag(rval)<-(1-prob)
    rval
}

multi.strat2Dcheck<-function(id,strata,probs){
   nstage<-ncol(id)
   rval<-vector("list",nstage)
   for(stage in 1:nstage){
       uid<-!duplicated(id[,stage])
       rval[[stage]]<-list(id=id[,stage],
                           dcheck=strat2Dcheck(strata[uid,stage], probs[uid,stage]
                           ))
   }
   rval
 }

overton2Dcheck<-function(prob,strat=rep(1,length(prob))){
    fbar<-outer(prob,prob,"+")/2
    n<-ave(strat,strat,FUN=length)
    rval<-(1-(n-fbar)/(n-1))
    rval[!outer(strat,strat,"==") | fbar==1]<-0
    diag(rval)<-(1-diag(fbar))
    rval
}

multi.overton2Dcheck<-function(id,strata,prob){
    nstage<-ncol(id)
    rval<-vector("list",nstage)
    for(stage in 1:nstage){
        uid<-!duplicated(id[,stage])
        rval[[stage]]<-list(id=id[,stage],
                            dcheck=overton2Dcheck(prob[uid,stage],
                              strata[uid,stage])
                           )
      }
    rval
}

htvar.list<-function(xcheck, Dcheck){
    rval<-sapply(Dcheck, function(stagei)
             {htvar.matrix(rowsum(xcheck,stagei$id),stagei$dcheck)})
    rval
}


htvar.matrix<-function(xcheck, Dcheck){
  if (is.null(dim(xcheck)))
    xcheck<-as.matrix(xcheck)
  apply(xcheck,2, function(xicheck)
        apply(xcheck,2, function(xjcheck)
              crossprod(xicheck, Dcheck%*%xjcheck)
              ))   
}

ygvar.matrix<-function(xcheck,Dcheck){
  ht<-htvar.matrix(xcheck,Dcheck)
  corr<-apply(xcheck,2, function(xicheck)
        apply(xcheck,2, function(xjcheck)
              sum(Dcheck%*%(xicheck*xjcheck))
              ))   
  ht-corr
}

twophaseDcheck<-function(Dcheck1,subset,Dcheck2){
  Dcheck1a<-Dcheck1[subset,subset]
  -Dcheck1a*Dcheck2+Dcheck1a+Dcheck2
}



