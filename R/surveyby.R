##
##  tables of statistics.
##

svyby<-function(formula, by, design, FUN,...,  keep.var=FALSE,keep.names=TRUE){

  if (inherits(by, "formula"))
    byfactors<-model.frame(by, design$variables, na.action=na.pass)
  else
    byfactors<-as.data.frame(by)
  
  byfactor<-do.call("interaction", byfactors)
  uniques <- which(!duplicated(byfactors))
  unwrap<-function(x) c(statistic=unclass(x),sd=sqrt(diag(as.matrix(attr(x,"var")))))
  
  if (keep.var)
    rval<-t(sapply(uniques, function(i) unwrap(FUN(formula,design[byfactor==byfactor[i],],...))))
  else {
    rval<-sapply(uniques, function(i) FUN(formula,design[byfactor==byfactor[i],],...))
    if (is.matrix(rval)) rval<-t(rval)
  }
  if (NCOL(rval)>1)
    rval<-cbind(byfactors[uniques,,drop=FALSE], rval)
  else
    rval <-cbind(byfactors[uniques,,drop=FALSE], statistic=rval)

  if (keep.names)
    rownames(rval)<-as.character(byfactor[uniques])

  rval<-rval[do.call("order",rval),]

  if (!keep.names)
    rownames(rval)<-1:NROW(rval)
  
  attr(rval,"call")<-sys.call()
  rval
}
