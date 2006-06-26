##
##  tables of statistics.
##

svyby<-function(formula, by, design, FUN,..., deff=FALSE, keep.var=TRUE,
                keep.names=TRUE,verbose=FALSE,vartype=c("se","cv","cvpct","var"),
                drop.empty.groups=TRUE, covmat=FALSE){

  if (inherits(by, "formula"))
    byfactors<-model.frame(by, design$variables, na.action=na.pass)
  else
    byfactors<-as.data.frame(by)

  if(covmat){
    if (!inherits(design,"svyrep.design"))
      stop("covmat=TRUE not implemented for this design type")
  }
  
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
        rval<-c(statistics=if(covmat) unclass(x[[1]]) else unclass(x))
        nvar<-length(rval)
        rval<-c(rval,c(se=SE(x),
                       cv=cv(x),
                       `cv%`=cv(x)*100,
                       var=SE(x)^2)[rep((nvartype-1)*(nvar),each=nvar)+(1:nvar)])
        if(!is.null(attr(x,"deff")))
          rval<-c(rval,DEff=deff(x))
        rval
      }

      ## In dire need of refactoring (or rewriting)
      ## but it seems to work.
      results<-lapply(uniques,
                      function(i){
                        if(verbose) print(as.character(byfactor[i]))
                        if (inherits(formula,"formula"))
                          data<-formula
                        else
                          data<-subset(formula, byfactor %in% byfactor[i])
                        if (covmat) {
                          FUN(data,
                              design[byfactor %in% byfactor[i],],
                              deff=deff,...,return.replicates=TRUE)
                        } else {
                          FUN(data,
                              design[byfactor %in% byfactor[i],],
                              deff=deff,...)
                        }
                      })
      rval<-t(sapply(results, unwrap))
      if (covmat) {
        replicates<-do.call(cbind,lapply(results,"[[","replicates"))
        covmat.mat<-svrVar(replicates,design$scale,design$rscales)
      }
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

  expand.index<-function(index,reps,x=FALSE){
    print(index)
    ns<-max(index)
    if (x){
      i<-matrix(1:(ns*reps),ncol=reps)
      rval<-t(i[index,])
      
    } else{
      i<-matrix(1:(ns*reps), ncol=reps, nrow=ns, byrow=TRUE)
      rval<- i[index,]
    }
    print(rval)
    as.vector(rval)
  }

  if(drop.empty.groups){
      if (keep.names)
          rownames(rval)<-paste(byfactor[uniques])
      rval<-rval[order(byfactor[uniques]),]
      if(covmat){
        i<-expand.index(order(byfactor[uniques]),nstats)
        covmat.mat<-covmat.mat[i,i]
      }
  } else {
      a<-do.call("expand.grid", lapply(byfactors,function(f) levels(as.factor(f))))
      a<-cbind(a,matrix(NA, ncol=nr, nrow=nrow(a)))
      names(a)<-names(rval)
      a[match(byfactor[uniques], levels(byfactor)),]<-rval
      rval<-a
      if (keep.names)
          rownames(rval)<-levels(byfactor)
      if(covmat){
        tmp<-matrix(ncol=nrow(a)*nstats,nrow=nrow(a)*nstats)
        i<-expand.index(match(byfactor[uniques], levels(byfactor)),nstats,TRUE)
        tmp[i,i]<-covmat.mat
        covmat.mat<-tmp
      }
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

  if(covmat)
    attr(rval,"var")<-covmat.mat
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

vcov.svyby<-function(object,...){
  rval<-attr(object,"var")
  if(is.null(rval)){
    warning("Only diagonal elements of vcov() available")
    se<-SE(object)
    if(length(se)>1)
      rval<-diag(se^2)
    else
      rval<-as.matrix(se^2)
  }
  rval
}
