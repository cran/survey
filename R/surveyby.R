##
##  tables of statistics.
##

svyby<-function(formula, by, design, FUN,..., deff=FALSE, keep.var=TRUE,
                keep.names=TRUE,verbose=FALSE,vartype=c("se","cv","cvpct","var"),
                drop.empty.groups=TRUE){

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
  
  if(missing(vartype)) vartype<-"se"
  vartype<-match.arg(vartype,several.ok=TRUE)
  nvartype<-pmatch(vartype,eval(formals(sys.function())$vartype))
  if(any(is.na(nvartype))) stop("invalid vartype")
  
  if (keep.var){
      unwrap <-function(x){
          rval<-c(statistics=unclass(x))
          nvar<-length(rval)
          rval<-c(rval,c(se=SE(x),
                         cv=cv(x),
                         `cv%`=cv(x)*100,
                         var=SE(x)^2)[rep((nvartype-1)*(nvar),each=nvar)+(1:nvar)])
          if(!is.null(attr(x,"deff")))
              rval<-c(rval,DEff=deff(x))
          rval
      }

                       
  
    rval<-t(sapply(uniques,
                   function(i) {
                       if(verbose) print(as.character(byfactor[i]))
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
                     if(verbose) print(as.character(byfactor[i]))
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

  nr<-NCOL(rval)
  nstats<-nr/(1+ keep.var*length(vartype) +deff)

              
  if (nr>1)
    rval<-cbind(byfactors[uniques,,drop=FALSE], rval)
  else
    rval <-cbind(byfactors[uniques,,drop=FALSE], statistic=rval)


  if(drop.empty.groups){
      if (keep.names)
          rownames(rval)<-paste(byfactor[uniques])
      rval<-rval[order(byfactor[uniques]),]
  } else {
      a<-do.call("expand.grid", lapply(byfactors,function(f) levels(as.factor(f))))
      a<-cbind(a,matrix(NA, ncol=nr, nrow=nrow(a)))
      names(a)<-names(rval)
      a[match(byfactor[uniques], levels(byfactor)),]<-rval
      rval<-a
      if (keep.names)
          rownames(rval)<-levels(byfactor)
  }
                  

  
  attr(rval,"svyby")<-list(margins=1:NCOL(byfactors),nstats=nstats,
                           vars=if(keep.var) length(vartype) else 0,
                           deffs=deff,
                           statistic=deparse(substitute(FUN)),
                           variables= names(rval)[-(1:NCOL(byfactors))][1:nstats],
                           vartype=vartype
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
    vartype<-attr(object,"svyby")$vartype
    if (pos<-match("se",vartype,0))
        object[,max(aa$margins)+aa$nstats*pos+(1:aa$nstats)]
    else if (pos<-match("var",vartype,0))
        sqrt(object[,max(aa$margins)+aa$nstats*pos+(1:aa$nstats)])
    else if (pos<-match("cv",vartype,0))
        object[,max(aa$margins)+aa$nstats*pos+(1:aa$nstats)]*coef(object)
    else if (pos<-match("cvpct",vartype,0))
         object[,max(aa$margins)+aa$nstats*pos+(1:aa$nstats)]*coef(object)/100
    else stop("This can't happen")
           
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
