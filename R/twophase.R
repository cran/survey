##
##
twophase<-function(id,strata=NULL, probs=NULL, weights=NULL, fpc=NULL,
                   subset, data){

  d1<-svydesign(id=id[[1]],strata=strata[[1]],weights=weights[[1]],
                probs=probs[[1]],fpc=fpc[[1]],data=data)

  if(inherits(subset,"formula"))
    subset<-eval.parent(model.frame(subset,data=data))[[1]]

  if(!is.logical(subset) && sort(unique(subset))==c(0,1))
      subset<-as.logical(subset)
  
  ## work out phase-two fpc
  if (is.null(fpc[[2]])){
    complete.vars<-names(data)[apply(data,2,function(v) all(!is.na(v)))]
    if (all(c(all.vars(id[[2]]),all.vars(strata[[2]])) %in% complete.vars)){
      dfpc<-svydesign(id=id[[2]],strata=strata[[2]],data=data,probs=NULL)
      popsize<-mapply(function(s,i) ave(!duplicated(i),s,FUN=sum), dfpc$strata, dfpc$cluster)
      rm(dfpc)
    } else {
      warning("Second-stage fpc not specified and not computable")
      popsize<-NULL
    }
  } else popsize<-NULL

  d2<-svydesign(id=id[[2]], strata=strata[[2]], probs=probs[[2]],
                weights=weights[[2]], fpc=fpc[[2]], data=data[subset,])

  ## ugly hack to get nicer labels
  if(!is.null(fpc[[2]])){
    d2call<-bquote(svydesign(id=.(id[[2]]),strata=.(strata[[2]]), probs=.(probs[[2]]),
                              weights=.(weights[[2]]), fpc=.(fpc[[2]])))
  } else{
    d2call<-bquote(svydesign(id=.(id[[2]]),strata=.(strata[[2]]), probs=.(probs[[2]]),
                              weights=.(weights[[2]]), fpc=`*phase1*`))
  }
  for(i in names(d2call)[-1])
    d2call[[i]]<-d2call[[i]]
  d2$call<-d2call
  d1call<-bquote(svydesign(id=.(id[[1]]),strata=.(strata[[1]]), probs=.(probs[[1]]),
                              weights=.(weights[[1]]), fpc=.(fpc[[1]])))
  for(i in names(d1call)[-1])
    d1call[[i]]<-d1call[[i]]
  d1$call<-d1call


  ## Add phase 2 fpc and probs if they were computed rather than specified.
  if (!is.null(popsize))
    d2$fpc<-as.fpc(popsize[subset,,drop=FALSE],d2$strata,d2$cluster)
  if(is.null(probs[[2]]) && is.null(weights[[2]]) && !is.null(d2$fpc$popsize)){
    d2$allprob<-1/weights(d2$fpc,final=FALSE)
    d2$prob<-apply(as.data.frame(d2$allprob),1,prod)
  }
  
  d2$variables<-NULL
  rval<-list(phase1=list(full=d1,sample=subset(d1,subset)),
             phase2=d2,
             subset=subset)
  rval$prob<-rval$phase1$sample$prob

  nunique<-function(x) sum(!duplicated(x))
  m<-NCOL(rval$phase1$sample$cluster)
  if(d2$has.strata){
      if (inherits(strata[[2]],"formula"))
          sa<-eval(attr(terms(strata[[2]]),"variables")[[2]],d1$variables)
      else
          sa<-d1$strata[,1]
      cm<-rval$phase1$full$cluster[,m]
      if (nunique(sa)!=nunique(sa[subset]))
          stop("Some phase-2 strata have zero sampling fraction")
      rval$usu<-ave(cm[subset],sa[subset],FUN=nunique)/ave(cm,sa,FUN=nunique)[subset]
  } else {
      rval$usu<-drop(with(rval$phase1$sample,ave(cluster[,m], strata[,m], FUN=nunique)/fpc$sampsize))
  }
  
  if (length(rval$phase1$sample$prob)==length(d2$prob))
    rval$prob<-rval$phase1$sample$prob*d2$prob
  else{
    rval$prob<-rep(Inf,length(rval$phase1$sample$prob))
    rval$prob[subset]<-rval$prob[subset]*d2$prob
  }
  rval$call<-sys.call()
  class(rval) <- c("twophase","survey.design")
  rval
}

print.twophase<-function(x,...){
  cat("Two-phase design:")
  print(x$call)
  cat("Phase 1:\n")
  print(x$phase1$full)
  cat("Phase 2:\n")
  print(x$phase2)
  invisible(x)
}

summary.twophase<-function(object,...){
  class(object)<-"summary.twophase"
  object
}

print.summary.twophase<-function(x,...,varnames=TRUE){
  cat("Two-phase design:")
  print(x$call)
   cat("Phase 1:\n")
  print(x$phase1$full,design.summaries=TRUE,varnames=FALSE)
  cat("Phase 2:\n")
  print(x$phase2,design.summaries=TRUE, varnames=FALSE)
  if (varnames){
    cat("Data variables:\n")
    print(names(x$phase1$full$variables))
  }
  invisible(x)
}

twophasevar<-function(x,design){
  d1 <- design$phase1$full
  u <- matrix(0,nrow=length(d1$prob),ncol=NCOL(x))
  u[design$subset,] <- x*sqrt(design$usu)
  u[is.na(u)]<-0
  vphase1 <- svyrecvar(u,d1$cluster, d1$strata, d1$fpc, postStrata=d1$postStrata)
  u2 <- x*sqrt(design$phase1$sample$prob) 
  if (NROW(u2)!=length(design$phase2$prob))
      u2 <- u2[subset,,drop=FALSE]
  else
      u2[is.na(u2)]<-0
  vphase2 <- with(design, svyrecvar(u2, phase2$cluster, phase2$strata,
                                  phase2$fpc, postStrata=phase2$postStrata))
  rval <- vphase1+vphase2
  attr(rval,"phases")<-list(phase1=vphase1, phase2=vphase2)
  rval
}


svytotal.twophase<-function(x,design, na.rm=FALSE, deff=FALSE,...){

  
  if (inherits(x,"formula")){
    ## do the right thing with factors
    mf<-model.frame(x,design$phase1$sample$variables,
                    na.action=na.pass)
    xx<-lapply(attr(terms(x),"variables")[-1],
               function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
    cols<-sapply(xx,NCOL)
    x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
    scols<-c(0,cumsum(cols))
    for(i in 1:length(xx)){
      x[,scols[i]+1:cols[i]]<-xx[[i]]
    }
    colnames(x)<-do.call("c",lapply(xx,colnames))
  } else if(typeof(x) %in% c("expression","symbol"))
    x<-eval(x, design$variables)
  
  x<-as.matrix(x)
  
  if (na.rm){
      nas<-rowSums(is.na(x))
      design<-design[nas==0,]
      if(length(nas)>length(design$prob))
          x<-x[nas==0,,drop=FALSE]
      else
          x[nas>0,]<-0
  }
  
  N<-sum(1/design$prob)
  total <- colSums(x/as.vector(design$prob),na.rm=na.rm)
  class(total)<-"svystat"
  attr(total, "var")<-v<-twophasevar(x/design$prob,design)
  attr(total,"statistic")<-"total"
  
  if (is.character(deff) || deff){
    nobs<-NROW(design$cluster)
    if (deff=="replace")
      vsrs<-svyvar(x,design,na.rm=na.rm)*sum(weights(design)^2)*(N-nobs)/N
    else
      vsrs<-svyvar(x,design,na.rm=na.rm)*sum(weights(design)^2)
    attr(total, "deff")<-v/vsrs
  }
  
  
  return(total)
}

svymean.twophase<-function(x,design, na.rm=FALSE,deff=FALSE,...){
  
  if (inherits(x,"formula")){
    ## do the right thing with factors
    mf<-model.frame(x,design$phase1$sample$variables
                    ,na.action=na.pass)
    xx<-lapply(attr(terms(x),"variables")[-1],
               function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
    cols<-sapply(xx,NCOL)
    x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
    scols<-c(0,cumsum(cols))
    for(i in 1:length(xx)){
      x[,scols[i]+1:cols[i]]<-xx[[i]]
    }
    colnames(x)<-do.call("c",lapply(xx,colnames))
  }
  else if(typeof(x) %in% c("expression","symbol"))
    x<-eval(x, design$phase1$sample)
  
  x<-as.matrix(x)
  
  if (na.rm){
    nas<-rowSums(is.na(x))
    design<-design[nas==0,]
    if(length(nas)>length(design$prob))
        x<-x[nas==0,,drop=FALSE]
    else
        x[nas>0,]<-0
  }
  
  pweights<-1/design$prob
  psum<-sum(pweights)
  average<-colSums(x*pweights/psum)
  x<-sweep(x,2,average)
  v<-twophasevar(x*pweights/psum,design)
  attr(average,"var")<-v
  attr(average,"statistic")<-"mean"
  class(average)<-"svystat"
  if (is.character(deff) || deff){
      nobs<-NROW(design$cluster)
      if(deff=="replace"){
        vsrs<-svyvar(x,design,na.rm=na.rm)/(nobs)
      } else {
        if(psum<nobs) {
          vsrs<-NA*v
          warning("Sample size greater than population size: are weights correctly scaled?")
        } else{
          vsrs<-svyvar(x,design,na.rm=na.rm)*(psum-nobs)/(psum*nobs)
        }
      }
      attr(average, "deff")<-v/vsrs
  }
  
  return(average)
}

model.frame.twophase<-function(formula,phase=2,...){
  if (phase==1)
    formula$phase1$full$variables
  else 
    formula$phase1$sample$variables
}

svyratio.twophase<-function(numerator, denominator, design, separate=FALSE,na.rm=FALSE,...){

    if (separate){
      strats<-sort(unique(design$phase2$strata[,1]))
      if (!design$phase2$has.strata)
        warning("Separate and combined ratio estimators are the same for unstratified designs")
      rval<-list(ratios=lapply(strats,
                   function(s) {
                     tmp<-svyratio(numerator, denominator,
                                   subset(design, design$phase2$strata[,1] %in% s),
                                   separate=FALSE,...)
                     tmp$call<-bquote(Stratum==.(s))
                     tmp}))
      names(rval$ratios)<-strats
   
      class(rval)<-c("svyratio_separate")
      rval$call<-sys.call()
      rval$strata<-strats
      return(rval)
    }
  
    if (inherits(numerator,"formula"))
        numerator<-model.frame(numerator,model.frame(design),na.action=na.pass)
    else if(typeof(numerator) %in% c("expression","symbol"))
        numerator<-eval(numerator, design$variables)
    if (inherits(denominator,"formula"))
        denominator<-model.frame(denominator,model.frame(design),na.action=na.pass)
    else if(typeof(denominator) %in% c("expression","symbol"))
        denominator<-eval(denominator, model.frame(design))

    nn<-NCOL(numerator)
    nd<-NCOL(denominator)

    all<-cbind(numerator,denominator)
    nas<-!complete.cases(all)
    if (na.rm){
        design<-design[!nas,]
        all<-all[!nas,,drop=FALSE]
        numerator<-numerator[!nas,,drop=FALSE]
        denominator<-denominator[!nas,,drop=FALSE]
    }
    allstats<-svytotal(all, design) 
    rval<-list(ratio=outer(allstats[1:nn],allstats[nn+1:nd],"/"))


    vars<-matrix(ncol=nd,nrow=nn)
    for(i in 1:nn){
      for(j in 1:nd){
        r<-(numerator[,i]-rval$ratio[i,j]*denominator[,j])/sum(denominator[,j]/design$prob)
        vars[i,j]<-twophase(r*1/design$prob, design)
      }
    }
    colnames(vars)<-names(denominator)
    rownames(vars)<-names(numerator)
    rval$var<-vars
    rval$call<-sys.call()
    class(rval)<-"svyratio"
    rval
    
  }


"[.twophase"<-function (x,i, ..., drop=TRUE){
  if (!missing(i)){ 
      if (is.calibrated(x$phase1$full) || is.calibrated(x$phase2) || !drop){
          ## Set weights to zero: no memory saving possible
          ## There should be an easier way to complement a subscript..
          if (is.logical(i)){
              x$prob[!i]<-Inf
              x$phase2$prob[!i]<-Inf
          } else if (is.numeric(i) && length(i)){
              x$prob[-i]<-Inf
              x$phase2$prob[-i]<-Inf
          } else {
              tmp<-x$prob[i,]
              x$prob<-rep(Inf, length(x$prob))
              x$prob[i,]<-tmp
          }
          index<-is.finite(x$prob)
          psu<-!duplicated(x$phase2$cluster[index,1])
          tt<-table(x$phase2$strata[index,1][psu])
          if(any(tt==1)){
              warning(sum(tt==1)," strata have only one PSU in this subset.")
          }
      } else {
          ## subset everything.
          x$prob<-x$prob[i]
          if (is.logical(i))
              x$subset[x$subset]<- i
          else if (is.numeric(i) && length(i))
              x$subset[which(x$subset)[-i]]<- FALSE
          else
              x$subset<-FALSE & x$subset
          x$usu<-x$usu[i]
          x$phase1$sample<-x$phase1$sample[i,...,drop=TRUE]
          x$phase2<-x$phase2[i,...,drop=TRUE]
        }
  } else {
    x$phase1$full<-x$phase1$full[,...]
    x$phase1$sample<-x$phase1$sample[,...]
    x$phase2<-x$phase2[,...]
  }
  x
}

dim.twophase<-function(x,...){
	dim(x$phase1$sample$variables)
}

na.fail.twophase<-function(object,...){
	tmp<-na.fail(object$phase1$sample$variables,...)
	object
}

na.omit.twophase<-function(object,...){
  tmp<-na.omit(object$phase1$sample$variables,...)
  omit<-attr(tmp,"na.action")
  if (length(omit)){
    object<-object[-omit,]
    object$phase1$sample$variables<-tmp
    attr(object,"na.action")<-omit
  }
  object
}

na.exclude.twophase<-function(object,...){
	tmp<-na.exclude(object$phase1$sample$variables,...)
	exclude<-attr(tmp,"na.action")
	if (length(exclude)){
           object<-object[-exclude,]
	   object$phase1$sample$variables<-tmp
	   attr(object,"na.action")<-exclude
	}
	object
}


update.twophase<-function(object,...){

  dots<-substitute(list(...))[-1]
  newnames<-names(dots)
  
  for(j in seq(along=dots)){
    object$phase1$sample$variables[,newnames[j]]<-eval(dots[[j]],object$phase1$sample$variables, parent.frame())
    object$phase1$full$variables[,newnames[j]]<-eval(dots[[j]],object$phase1$full$variables, parent.frame())
  }
  
  object$call<-sys.call(-1)
  object 
}

subset.twophase<-function(x,subset,...){
        e <- substitute(subset)
        r <- eval(e, x$phase1$sample$variables, parent.frame())
        r <- r & !is.na(r) 
        x<-x[r,]
	x$call<-sys.call(-1)
	x
}


calibrate.twophase<-function(design, phase=2, formula, population,...){

    if (phase==1){
        
        phase1<-calibrate(design$phase1$full,formula, population, ...)
        design$phase1$full<-phase1
        design$phase1$sample<-phase1[design$subset,]
        
    } else if(phase==2){

        if (missing(population) || is.null(population)){
            ## calibrate to phase 1 totals
            population<-colSums(model.matrix(formula,
                               model.frame(formula, design$phase1$full$variables)))
        }
        
        phase2<-design$phase2
        phase2$variables<-design$phase1$sample$variables
        phase2<-calibrate(phase2,formula,population,...)
        g<-design$phase2$prob/phase2$prob
        phase2$variables<-NULL
        design$phase2<-phase2
	design$usu<-design$usu/g

    } else stop("`phase' must be 1 or 2")

    
    if (length(design$phase1$sample$prob)==length(design$phase2$prob))
        design$prob<-design$phase1$sample$prob*design$phase2$prob
    else{
        design$prob<-rep(Inf,length(design$phase1$sample$prob))
        design$prob[subset]<-design$prob[subset]*design$phase2$prob
    }

    design$call<-sys.call(-1)
    
    design

}


postStratify.twophase<-function(design, ...) {
	stop("postStratify not yet implemented for two-phase designs. Use calibrate()")
}

