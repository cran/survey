##
##  tables of statistics.
##
svyby<-function(formula, by, design,...) UseMethod("svyby",design)


subset_drops_rows<-function(design) !(is.pps(design) || is.calibrated(design))

unwrap<-function(x,nvartype) UseMethod("unwrap",x)
unwrap2<-function(x) UseMethod("unwrap2",x)

unwrap2.default <- function(x){
    if(!is.null(attr(x, "deff")))
        c(statistic = unclass(x),
          DEff = deff(x))
    else c(statistic = unclass(x))
}


unwrap2.newsvyquantile <- function(x){
    rval<-do.call(rbind,x)
    rownames(rval)<-names(x)
    rval
}



unwrap.default <-function(x,nvartype){
    rval<-c(coef(x))
    nvar<-length(rval)
    rval<-c(rval,c(se=SE(x),
                   ci_l=confint(x)[,1],
                   ci_u=confint(x)[,2],
                   cv=cv(x,warn=FALSE),
                   `cv%`=cv(x,warn=FALSE)*100,
                   var=SE(x)^2)[rep((nvartype-1)*(nvar),each=nvar)+(1:nvar)])
    if(!is.null(attr(x,"deff")))
        rval<-c(rval,DEff=deff(x))
    rval
}




strings_to_factors<-function(formula,  design){
    allv<-intersect(all.vars(formula), colnames(design))
    vclass<-sapply(model.frame(design)[,allv,drop=FALSE], class)
    if (!any(vclass=="character")) return(design)
    vfix<-names(vclass)[vclass=="character"]
    l<-as.list(vfix)
    names(l)<-vfix
    fl<-lapply(l, function(li) bquote(factor(.(as.name(li)))))
    expr<-bquote(update(design, ..(fl)), splice=TRUE)
    eval(expr)
}


svyby.default<-function(formula, by, design, FUN,..., deff=FALSE, keep.var=TRUE,
                        keep.names=TRUE, verbose=FALSE, vartype=c("se","ci","ci","cv","cvpct","var"),
                        drop.empty.groups=TRUE, covmat=FALSE, return.replicates=FALSE, na.rm.by=FALSE,
                        na.rm.all=FALSE, stringsAsFactors=TRUE,
                        multicore=getOption("survey.multicore")){

  if (inherits(by, "formula"))
    byfactors<-model.frame(by, model.frame(design), na.action=na.pass)
  else
    byfactors<-as.data.frame(by)

  if(return.replicates){
    if (!inherits(design,"svyrep.design"))
      stop("return.replicates=TRUE not implemented for this design type")
  }
    
    if(stringsAsFactors){
        design<-strings_to_factors(formula,design)
    }
    
  if (multicore && !requireNamespace("parallel",quietly=TRUE))
    multicore<-FALSE

  ## some people insist on using vectors rather than formulas
  ## so I suppose we should be nice to them
  if (!inherits(formula, "formula")){
      if (NROW(formula)!=NROW(byfactors))
          stop("'formula' is the wrong length")
      if (!(is.data.frame(formula) ||
            is.matrix(formula) ||
            is.vector(formula))){
          stop("invalid type for 'formula'")
      }
  }

  hasdeff<- is.character(deff) || deff
  
  ## all combinations that actually occur in this design
  byfactor<-do.call("interaction", byfactors)
  dropped<- weights(design,"sampling")==0
  if (na.rm.by) dropped<-dropped | apply(byfactors, 1, function(x) any(is.na(x)))
  if (na.rm.all){
    if (inherits(formula,"formula"))
      allx<-model.frame(formula,model.frame(design),na.action=na.pass)
    else
      allx<-formula
    dropped <- dropped | (!complete.cases(allx))
  }
  uniquelevels<-sort(unique(byfactor[!dropped]))
  uniques <- match(uniquelevels, byfactor)

  
  if(missing(vartype)) vartype<-"se"
  vartype<-match.arg(vartype,several.ok=TRUE)
  nvartype<-base::which(eval(formals(sys.function())$vartype) %in% vartype)
  if(any(is.na(nvartype))) stop("invalid vartype")
  
  if (keep.var){

      ## In dire need of refactoring (or rewriting)
      ## but it seems to work.
      results<-(if (multicore) parallel::mclapply else lapply)(uniques,
                      function(i){
                        if(verbose && !multicore) print(as.character(byfactor[i]))
                        if (inherits(formula,"formula"))
                          data<-formula
                        else
                          data<-subset(formula, byfactor %in% byfactor[i])
                        if (covmat || return.replicates) {
                          FUN(data,
                              design[byfactor %in% byfactor[i],],
                              deff=deff,...,return.replicates=TRUE)
                        } else {
                          FUN(data,
                              design[byfactor %in% byfactor[i],],
                              deff=deff,...)
                        }
                      })
      rval<-t(sapply(results, unwrap,nvartype=nvartype))
      if ((covmat && inherits(design, "svyrep.design")) || return.replicates) {
          replicates<-do.call(cbind,lapply(results,"[[","replicates"))
          attr(replicates,"scale")<-design$scale
          attr(replicates, "rscales")<-design$rscales
          attr(replicates, "mse")<-design$mse
        colnames(replicates)<-rep(as.character(uniquelevels), each=NCOL(replicates)/length(uniquelevels))
        covmat.mat<-svrVar(replicates,design$scale,design$rscales, mse=design$mse,coef=as.vector(sapply(results,coef)))
      } else{
        covmats<-lapply(results,vcov)
        ncovmat<-sum(sapply(covmats,ncol))
        covmat.mat<-matrix(0,ncol=ncovmat,nrow=ncovmat)
        j<-0
        for(i in 1:length(covmats)){
          ni<-nrow(covmats[[i]])
          covmat.mat[j+(1:ni),j+(1:ni)]<-covmats[[i]]
          j<-j+ni
        }
      }      
    } else {
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
  nstats<-nr/(1+ keep.var*(length(vartype)+ ("ci" %in% vartype)) + hasdeff)

              
  if (nr>1)
    rval<-cbind(byfactors[uniques,,drop=FALSE], rval)
  else
    rval <-cbind(byfactors[uniques,,drop=FALSE], statistic=rval)

  expand.index<-function(index,reps,x=FALSE){
    ns<-max(index)
    if (x){
      i<-matrix(1:(ns*reps),ncol=reps)
      rval<-t(i[index,])
      
    } else{
      i<-matrix(1:(ns*reps), ncol=reps, nrow=ns, byrow=TRUE)
      rval<- i[index,]
    }
    as.vector(rval)
  }

  if(drop.empty.groups){
      if (keep.names)
          rownames(rval)<-paste(byfactor[uniques])
      rval<-rval[order(byfactor[uniques]),]

      i<-expand.index(order(byfactor[uniques]),nstats)
      if (keep.var)
        covmat.mat<-covmat.mat[i,i]

  } else {
      a<-do.call("expand.grid", lapply(byfactors,function(f) levels(as.factor(f))))
      a<-cbind(a,matrix(NA, ncol=nr, nrow=nrow(a)))
      names(a)<-names(rval)
      a[match(byfactor[uniques], levels(byfactor)),]<-rval
      rval<-a
      if (keep.names)
          rownames(rval)<-levels(byfactor)
      if (keep.var){
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
  if (return.replicates)
    attr(rval,"replicates")<-replicates
  attr(rval,"call")<-sys.call()
  class(rval)<-c("svyby","data.frame")
  rval
}

SE.svyby <-function(object,...){
    aa<-attr(object,"svyby")
    if (!aa$vars) stop("Object does not contain variances")
    vartype<-attr(object,"svyby")$vartype
    vartype <- c("se","ci","ci","cv","cvpct","var")[c("se","ci","ci","cv","cvpct","var") %in% vartype]
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

coef.svyby<-function (object, ...)
{
    aa <- attr(object, "svyby")
    rval <- object[, max(aa$margins) + (1:aa$nstats)]
    if (is.null(dim(rval))){
       names(rval) <- row.names(object)
    } else {
        rval<-as.vector(as.matrix(rval))
        names(rval)<-outer(rownames(object),
                           gsub("statistics\\.","",aa$variables), paste, sep=":")
    }
    rval
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
    if (is.data.frame(se)) se<-as.vector(as.matrix(se))
    if(length(se)>1)
      rval<-diag(se^2)
    else
      rval<-as.matrix(se^2)
  }
  nms<-names(coef(object))
  dimnames(rval)<-list(nms,nms)
  rval
}


## allows for influence functions. Does not allow for replicates
	
svyby.survey.design2<-function(formula, by, design, FUN,..., deff=FALSE, keep.var=TRUE,
                        keep.names=TRUE, verbose=FALSE, vartype=c("se","ci","ci","cv","cvpct","var"),
                        drop.empty.groups=TRUE, covmat=FALSE, influence=covmat, na.rm.by=FALSE,
                        na.rm.all=FALSE, stringsAsFactors=TRUE,
                        multicore=getOption("survey.multicore")){

  if (inherits(by, "formula"))
    byfactors<-model.frame(by, model.frame(design), na.action=na.pass)
  else
    byfactors<-as.data.frame(by)

    if (stringsAsFactors){
        design<-strings_to_factors(formula, design)
    }

  if (multicore && !requireNamespace("parallel",quietly=TRUE))
    multicore<-FALSE

  ## some people insist on using vectors rather than formulas
  ## so I suppose we should be nice to them
  if (!inherits(formula, "formula")){
      if (NROW(formula)!=NROW(byfactors))
          stop("'formula' is the wrong length")
      if (!(is.data.frame(formula) ||
            is.matrix(formula) ||
            is.vector(formula))){
          stop("invalid type for 'formula'")
      }
  }

  hasdeff<- is.character(deff) || deff
  
  ## all combinations that actually occur in this design
  byfactor<-do.call("interaction", byfactors)
  dropped<- weights(design,"sampling")==0
  if (na.rm.by) dropped<-dropped | apply(byfactors, 1, function(x) any(is.na(x)))
  if (na.rm.all){
    if (inherits(formula,"formula"))
      allx<-model.frame(formula,model.frame(design),na.action=na.pass)
    else
      allx<-formula
    dropped <- dropped | (!complete.cases(allx))
  }
  uniquelevels<-sort(unique(byfactor[!dropped]))
  uniques <- match(uniquelevels, byfactor)

  
  if(missing(vartype)) vartype<-"se"
  vartype<-match.arg(vartype,several.ok=TRUE)
  nvartype<-base::which(eval(formals(sys.function())$vartype) %in% vartype)
  if(any(is.na(nvartype))) stop("invalid vartype")
  
  if (keep.var){
      ## In dire need of refactoring (or rewriting)
      ## but it seems to work.
      results<-(if (multicore) parallel::mclapply else lapply)(uniques,
          function(i){
              idx<- byfactor %in% byfactor[i]  ##save this for assembling influence functions
              if(verbose && !multicore) print(as.character(byfactor[i]))
              if (inherits(formula,"formula"))
                  data<-formula
              else
                  data<-subset(formula, byfactor %in% byfactor[i]) 
              if (covmat||influence ) {
                  r<-FUN(data,
                      design[byfactor %in% byfactor[i],],
                      deff=deff,...,influence=influence)
              } else {
                  r<-FUN(data,
                      design[byfactor %in% byfactor[i],],
                      deff=deff,...)
              }
              attr(r,"index")<-idx
              r
          })
      rval<-t(sapply(results, unwrap,nvartype=nvartype))
      if (covmat || influence) {
          ## do the influence function thing here
          ## have to handle both subset-> zero wt and subset -> gone
          infs<-lapply(results,attr, "influence")
          idxs<-lapply(results,attr, "index")
          if (all(sapply(infs,is.null)))
              stop("FUN does not return influence functions")
          inflmats<-vector("list",length(infs))
          for(i in seq_along(infs)){
              inflmats[[i]]<-matrix(0, ncol=NCOL(infs[[i]]),nrow=length(idxs[[i]]))
              if (subset_drops_rows(design)){
                  inflmats[[i]][idxs[[i]],]<-infs[[i]]
              } else {
                   inflmats[[i]][idxs[[i]],]<-infs[[i]][idxs[[i]],]
              }
          }
          inflmat<-do.call(cbind,inflmats)
          covmat.mat<-svyrecvar(inflmat,design$cluster,
                                     design$strata, design$fpc,
                                   postStrata=design$postStrata)
       } else{
        covmats<-lapply(results,vcov)
        ncovmat<-sum(sapply(covmats,ncol))
        covmat.mat<-matrix(0,ncol=ncovmat,nrow=ncovmat)
        j<-0
        for(i in 1:length(covmats)){
          ni<-nrow(covmats[[i]])
          covmat.mat[j+(1:ni),j+(1:ni)]<-covmats[[i]]
          j<-j+ni
        }
      }      
    } else {

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
  nstats<-nr/(1+ keep.var*(length(vartype)+ ("ci" %in% vartype)) + hasdeff)

              
  if (nr>1)
    rval<-cbind(byfactors[uniques,,drop=FALSE], rval)
  else
    rval <-cbind(byfactors[uniques,,drop=FALSE], statistic=rval)

  expand.index<-function(index,reps,x=FALSE){
    ns<-max(index)
    if (x){
      i<-matrix(1:(ns*reps),ncol=reps)
      rval<-t(i[index,])
      
    } else{
      i<-matrix(1:(ns*reps), ncol=reps, nrow=ns, byrow=TRUE)
      rval<- i[index,]
    }
    as.vector(rval)
  }

  if(drop.empty.groups){
      if (keep.names)
          rownames(rval)<-paste(byfactor[uniques])
      rval<-rval[order(byfactor[uniques]),]

      i<-expand.index(order(byfactor[uniques]),nstats)
      if (keep.var)
        covmat.mat<-covmat.mat[i,i]

  } else {
      a<-do.call("expand.grid", lapply(byfactors,function(f) levels(as.factor(f))))
      a<-cbind(a,matrix(NA, ncol=nr, nrow=nrow(a)))
      names(a)<-names(rval)
      a[match(byfactor[uniques], levels(byfactor)),]<-rval
      rval<-a
      if (keep.names)
          rownames(rval)<-levels(byfactor)
      if (keep.var){
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
  if (influence)
    attr(rval,"influence")<-inflmat
  attr(rval,"call")<-sys.call()
  class(rval)<-c("svyby","data.frame")
  rval
}


svybys<-function(formula,  bys,  design, FUN, ...){
  tms <- attr(terms(bys),"variables")[-1]	
  
  lapply(tms, function(tm){
      eval(bquote(svyby(.(formula),by=~.(tm),
                        design=.(design), FUN=.(FUN), ...)))
  })
  
}
