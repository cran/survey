##
##  tables of statistics.
##

svyby<-function(formula, by, design, FUN,..., deff=FALSE, keep.var=FALSE,
                keep.names=TRUE){

  if (inherits(by, "formula"))
    byfactors<-model.frame(by, design$variables, na.action=na.pass)
  else
    byfactors<-as.data.frame(by)
  
  byfactor<-do.call("interaction", byfactors)
  uniques <- which(!duplicated(byfactors))

  
  if (keep.var){
      unwrap <- function(x){
          if(!is.null(attr(x, "deff")))
              c(statistic = unclass(x),
                SE = sqrt(diag(as.matrix(attr(x, "var")))),
                DEff = deff(x))
          else c(statistic = unclass(x),
                 SE = sqrt(diag(as.matrix(attr(x, "var")))))
      }

    rval<-t(sapply(uniques,
                   function(i) unwrap(FUN(formula,design[byfactor %in% byfactor[i],],deff=deff,...))))
  } else {
      unwrap2 <- function(x){
          if(!is.null(attr(x, "deff")))
              c(statistic = unclass(x),
                DEff = deff(x))
          else c(statistic = unclass(x))
      }
      rval<-sapply(uniques,
                   function(i) unwrap2(FUN(formula,design[byfactor %in% byfactor[i],],deff=deff,...)))
      if (is.matrix(rval)) rval<-t(rval)
  }

  nstats<-NCOL(rval)/(1+keep.var+deff)
  
  if (NCOL(rval)>1)
    rval<-cbind(byfactors[uniques,,drop=FALSE], rval)
  else
    rval <-cbind(byfactors[uniques,,drop=FALSE], statistic=rval)

  if (keep.names)
    rownames(rval)<-paste(byfactor[uniques])

  rval<-rval[do.call("order",rval),]

  attr(rval,"svyby")<-list(margins=1:NCOL(byfactors),nstats=nstats,
                           vars=keep.var,
                           deffs=deff,
                           statistic=deparse(substitute(FUN)),
                           variables= names(rval)[-(1:NCOL(byfactors))][1:nstats]
                           )
  if (!keep.names)
    rownames(rval)<-1:NROW(rval)
  
  attr(rval,"call")<-sys.call()
  class(rval)<-c("svyby","data.frame")
  rval
}