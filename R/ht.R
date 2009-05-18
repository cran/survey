##not used yet
pi2Dcheck<-function(pmat,tolerance){
    rval<-(pmat-outer(diag(pmat)))/pmat
    rval[abs(rval)<tolerance]<-0
    as(rval,"sparseMatrix")
}

iidDcheck<-function(n){
  diag(n)
}

## not used yet: Overton's approximation for PPS
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

## used in twophase2var()
htvar.matrix<-function(xcheck, Dcheck){
  if (is.null(dim(xcheck)))
    xcheck<-as.matrix(xcheck)
  rval<-apply(xcheck,2, function(xicheck)
              apply(xcheck,2, function(xjcheck)
                    as.matrix(Matrix::crossprod(xicheck, Dcheck%*%xjcheck))
                    ))
  if(is.null(dim(rval))) dim(rval)<-c(1,1)
  rval
}

## not yet used.
ygvar.matrix<-function(xcheck,Dcheck){
  ht<-htvar.matrix(xcheck,Dcheck)
  if (is.null(dim(xcheck))){
    corr <- sum(Dcheck%*%(xcheck*xcheck))    
  } else {
    corr <- apply(xcheck,2, function(xicheck)
                apply(xcheck,2, function(xjcheck)
                      sum(Dcheck%*%(xicheck*xjcheck))
                      ))
  }
  ht-corr
}


