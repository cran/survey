##
##  tables of statistics.
##

svyby<-function(formula, by, design, FUN,..., deff=FALSE, keep.var=TRUE,
                keep.names=TRUE){

  if (inherits(by, "formula"))
    byfactors<-model.frame(by, design$variables, na.action=na.pass)
  else
    byfactors<-as.data.frame(by)
  
  byfactor<-do.call("interaction", byfactors)
  uniques <- which(!duplicated(byfactors))

  ## some people insist on using vectors rather than formulas
  ## so I suppose we should be nice to them
  if (!inherits(formula, "formula")){
      if (NROW(formula)!=length(byfactor))
          stop("'formula' is the wrong length")
      if (!(is.data.frame(formula) ||
            is.matrix(formula) ||
            is.vector(formula))){
          stop("invalid type for 'formula'")
      }
      
  }
  
  
  if (keep.var){
      unwrap <- function(x){
          if(!is.null(attr(x, "deff")))
              c(statistic = unclass(x),
                SE = SE(x),
                DEff = deff(x))
          else c(statistic = unclass(x),
                 SE = SE(x))
      }

    rval<-t(sapply(uniques,
                   function(i) {
                       if (inherits(formula,"formula"))
                           data<-formula
                       else
                           data<-subset(formula, byfactor %in% byfactor[i])
                       unwrap(FUN(data,
                                  design[byfactor %in% byfactor[i],],
                                  deff=deff,...)) }
                   ))
  } else {
      unwrap2 <- function(x){
          if(!is.null(attr(x, "deff")))
              c(statistic = unclass(x),
                DEff = deff(x))
          else c(statistic = unclass(x))
      }
      rval<-sapply(uniques,
                   function(i) {
                       if (inherits(formula,"formula"))
                           data<-formula
                       else
                           data<-subset(formula, byfactor %in% byfactor[i])
                       unwrap2(FUN(data,
                                   design[byfactor %in% byfactor[i],],
                                   deff=deff,...))}
                   )
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

SE.svyby <-function(object,...){
    aa<-attr(object,"svyby")
    if (!aa$vars) stop("Object does not contain variances")
    object[,max(aa$margins)+aa$nstats+(1:aa$nstats)]
}

coef.svyby<-function(object,...){
    aa<-attr(object,"svyby")
    object[,max(aa$margins)+(1:aa$nstats)]
}

deff.svyby<-function(object,...){
    aa<-attr(object,"svyby")
    if (!aa$deffs) stop("object does not have design effect information")
    object[,max(aa$margins)+aa$nstats*(1+aa$vars)+(1:aa$nstats)]
}
