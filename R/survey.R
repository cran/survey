svydesign<-function(ids,probs=NULL,strata=NULL,variables=NULL, fpc=NULL,
                    data=NULL, nest=FALSE, check.strata=!nest,weights=NULL){

    ## less memory-hungry version for sparse tables
    interaction<-function (..., drop = TRUE) {
        args <- list(...)
        narg <- length(args)
        if (narg == 1 && is.list(args[[1]])) {
            args <- args[[1]]
            narg <- length(args)
        }
        
        ls<-sapply(args,function(a) length(levels(a)))
        ans<-do.call("paste",c(lapply(args,as.character),sep="."))
        ans<-factor(ans)
        return(ans)
        
    }


    na.failsafe<-function(object,...){
      if (NCOL(object)==0)
        object
      else na.fail(object)
    }
    
     if(inherits(ids,"formula")) {
	 mf<-substitute(model.frame(ids,data=data,na.action=na.failsafe))   
	 ids<-eval.parent(mf)
	} else if (!is.null(ids))
            ids<-na.fail(data.frame(ids))

     if(inherits(probs,"formula")){
	mf<-substitute(model.frame(probs,data=data,na.action=na.failsafe))
	probs<-eval.parent(mf)
	}
     
     if(inherits(weights,"formula")){
       mf<-substitute(model.frame(weights,data=data,na.action=na.failsafe))
       weights<-eval.parent(mf)
     } else if (!is.null(weights))
         weights<-na.fail(data.frame(weights))
     
     if(!is.null(weights)){
       if (!is.null(probs))
         stop("Can't specify both sampling weights and probabilities")
       else
         probs<-1/weights
     }

      

    if (!is.null(strata)){
      if(inherits(strata,"formula")){
        mf<-substitute(model.frame(strata,data=data, na.action=na.failsafe))
        strata<-eval.parent(mf)
      }
      if(is.list(strata))
        strata<-na.fail(do.call("interaction", strata))
      if (!is.factor(strata))
        strata<-factor(strata)
      has.strata<-TRUE
    } else {
      strata<-factor(rep(1,NROW(ids)))
      has.strata <-FALSE
    }
    
    if (inherits(variables,"formula")){
        mf<-substitute(model.frame(variables,data=data))
        variables <- eval.parent(mf)
    } else if (is.null(variables)){
        variables<-data
    } else
        variables<-data.frame(variables)

    
     if (inherits(fpc,"formula")){
       mf<-substitute(model.frame(fpc,data=data,na.action=na.failsafe))
       fpc<-eval.parent(mf)
       if (length(fpc))
         fpc<-fpc[,1]
     }
     
    if (is.null(ids) || NCOL(ids)==0)
	ids<-data.frame(.id=seq(length=NROW(variables)))

     ## force subclusters nested in clusters
     if (nest && NCOL(ids)>1){
      N<-ncol(ids)
      for(i in 2:(N)){
          ids[,i]<-do.call("interaction", ids[,1:i,drop=TRUE])
      }
    }
     ## force clusters nested in strata
     if (nest && has.strata && NCOL(ids)){
       N<-NCOL(ids)
       for(i in 1:N)
         ids[,i]<-do.call("interaction", list(strata, ids[,i]))
     }

    ## check if clusters nested in strata 
     if (check.strata && nest)
      warning("No point in check.strata=TRUE if nest=TRUE")
    if(check.strata && !is.null(strata) && NCOL(ids)){
       sc<-rowSums(table(ids[,1],strata)>0)
       if(any(sc>1)) stop("Clusters not nested in strata")
    }

    ## Put degrees of freedom (# of PSUs in each stratum) in object, to 
    ## allow subpopulations
    if (NCOL(ids)){
        nPSU<-table(strata[!duplicated(ids[,1])])
    }


     if (!is.null(fpc)){

       if (NCOL(ids)>1){
         if (all(fpc<1))
           warning("FPC is not currently supported for multi-stage sampling")
         else
           stop("Can't compute FPC from population size for multi-stage sampling")
       }
       
       ## Finite population correction: specified per observation
       if (is.numeric(fpc) && length(fpc)==NROW(variables)){
         tbl<-by(fpc,list(strata),unique)
         if (any(sapply(tbl,length)!=1))
           stop("fpc not constant within strata")
         fpc<-data.frame(strata=factor(rownames(tbl),levels=levels(strata)),
                         N=as.vector(tbl))
       }
       ## Now reduced to fpc per stratum
       nstr<-table(strata[!duplicated(ids[[1]])])
       
       if (all(fpc[,2]<=1)){
         fpc[,2]<- nstr[match(as.character(fpc[,1]), names(nstr))]/fpc[,2]
       } else if (any(fpc[,2]<nstr[match(as.character(fpc[,1]), names(nstr))]))
         stop("Over 100% sampling in some strata")
       
     }

    ## if FPC specified, but no weights, use it for weights
    if (is.null(probs) && is.null(weights) && !is.null(fpc)){
      pstr<-nstr[match(as.character(fpc[,1]), names(nstr))]/fpc[,2]
      probs<-pstr[match(as.character(strata),as.character(fpc[,1]))]
      probs<-as.vector(probs)
    }

    
    if (is.numeric(probs) && length(probs)==1)
        probs<-rep(probs, NROW(variables))
    
    if (length(probs)==0) probs<-rep(1,NROW(variables))
    
    if (NCOL(probs)==1) probs<-data.frame(probs)

    rval<-list(cluster=ids)
    rval$strata<-strata
    rval$has.strata<-has.strata
    rval$prob<- apply(probs,1,prod) 
    rval$allprob<-probs
    rval$call<-match.call()
    rval$variables<-variables
    rval$fpc<-fpc
    rval$call<-sys.call()
    rval$nPSU<-nPSU
    class(rval)<-"survey.design"
    rval
  }

print.survey.design<-function(x,varnames=FALSE,design.summaries=FALSE,...){
  n<-NROW(x$cluster)
  if (x$has.strata) cat("Stratified ")
  un<-length(unique(x$cluster[,1]))
  if(n==un){
    cat("Independent Sampling design\n")
    is.independent<-TRUE
  } else {
    cat(NCOL(x$cluster),"- level Cluster Sampling design\n")
    nn<-lapply(x$cluster,function(i) length(unique(i)))
    cat(paste("With (",paste(unlist(nn),collapse=","),") clusters.\n"))
    is.independent<-FALSE
  }
  print(x$call)
  if (design.summaries){
    cat("Probabilities:\n")
    print(summary(x$prob))
    if(x$has.strata){
      cat("Stratum sizes: \n")
      a<-rbind(obs=table(x$strata),
	       design.PSU=x$nPSU,
               actual.PSU=if(!is.independent || !is.null(x$fpc))
               table(x$strata[!duplicated(x$cluster[,1])]))
      print(a)
    }
    if (!is.null(x$fpc)){
      if (x$has.strata) {
        cat("Population stratum sizes (PSUs): \n")
        print(x$fpc)
      } else {
        cat("Population size (PSUs):",x$fpc[,2],"\n")
      }
    }
  }
  if (varnames){
    cat("Data variables:\n")
    print(names(x$variables))
  }
  invisible(x)
}

"[.survey.design"<-function (x,i, ...){
  
  if (!missing(i)){ 
    x$variables<-"[.data.frame"(x$variables,i,...,drop=FALSE)
    x$cluster<-x$cluster[i,,drop=FALSE]
    x$prob<-x$prob[i]
    x$allprob<-x$allprob[i,,drop=FALSE]
    x$strata<-x$strata[i]
  } else {
    x$variables<-x$variables[,...]
  }
  
  x
}

"[<-.survey.design"<-function(x, ...,value){
  if (inherits(value, "survey.design"))
    value<-value$variables
  x$variables[...]<-value
  x
}

dim.survey.design<-function(x,...){
	dim(x$variables)
}

na.fail.survey.design<-function(object,...){
	tmp<-na.fail(object$variables,...)
	object
}

na.omit.survey.design<-function(object,...){
  tmp<-na.omit(object$variables,...)
  omit<-attr(tmp,"na.action")
  if (length(omit)){
    object$cluster<-object$cluster[-omit,,drop=FALSE]
    object$prob<-object$prob[-omit]
    object$allprob<-object$allprob[-omit,,drop=FALSE]
    object$variables<-tmp
    attr(object,"na.action")<-omit
  }
  object
}

na.exclude.survey.design<-function(object,...){
	tmp<-na.exclude(object$variables,...)
	exclude<-attr(tmp,"na.action")
	if (length(exclude)){
	   object$cluster<-object$cluster[-exclude,,drop=FALSE]
	   object$prob<-object$prob[-exclude]
	   object$allprob<-object$allprob[-exclude,,drop=FALSE]
	   object$variables<-tmp
	   attr(object,"na.action")<-exclude
	}
	object
}


update.survey.design<-function(object,...){

  dots<-substitute(list(...))[-1]
  newnames<-names(dots)
  
  for(j in seq(along=dots)){
    object$variables[,newnames[j]]<-eval(dots[[j]],object$variables, parent.frame())
  }
  
  object$call<-sys.call()
  object 
}

subset.survey.design<-function(x,subset,...){
        e <- substitute(subset)
        r <- eval(e, x$variables, parent.frame())
        r <- r & !is.na(r) 
        x<-x[r,]
	x$call<-sys.call()
	x
}

summary.survey.design<-function(object,...){
  class(object)<-"summary.survey.design"
  object
}

print.summary.survey.design<-function(x,...){
  y<-x
  class(y)<-"survey.design"
  print(y,varnames=TRUE,design.summaries=TRUE,...)
}	
     
svyCprod<-function(x, strata, psu, fpc, nPSU,
                   lonely.psu=getOption("survey.lonely.psu")){

  x<-as.matrix(x)
  n<-NROW(x)

  ##First collapse over PSUs

  if (is.null(strata)){
    strata<-rep("1",n)
    if (!is.null(nPSU))
        names(nPSU)<-"1"
  }
  else
    strata<-as.character(strata) ##can't use factors as indices in for()'

  if (!is.null(psu)){
    x<-rowsum(x, psu, reorder=FALSE)
    strata<-strata[!duplicated(psu)]
    n<-NROW(x)
  }
  
  if (!is.null(nPSU)){
      obsn<-table(strata)
      dropped<-nPSU[match(names(obsn),names(nPSU))]-obsn
      if(sum(dropped)){
        xtra<-matrix(0,ncol=NCOL(x),nrow=sum(dropped))
        strata<-c(strata,rep(names(dropped),dropped))
      	if(is.matrix(x))
	   x<-rbind(x,xtra)
        else
	   x<-c(x,xtra)
        n<-NROW(x)
      }
  }

  if(is.null(strata)){
      x<-t(t(x)-colMeans(x))
  } else {
      strata.means<-drop(rowsum(x,strata, reorder=FALSE))/drop(rowsum(rep(1,n),strata, reorder=FALSE))
      if (!is.matrix(strata.means))
          strata.means<-matrix(strata.means, ncol=NCOL(x))
      x<- x- strata.means[ match(strata, unique(strata)),,drop=FALSE]
  }
  
  p<-NCOL(x)
  v<-matrix(0,p,p)
  
  ss<-unique(strata)
  for(s in ss){
      this.stratum <- strata %in% s
      
      ## original number of PSUs in this stratum 
      ## before missing data/subsetting
      this.n <-nPSU[match(s,names(nPSU))]
      
      this.df <- this.n/(this.n-1)	
      
      if (is.null(fpc))
          this.fpc <- 1
      else{
          this.fpc <- fpc[,2][ fpc[,1]==as.character(s)]
          this.fpc <- (this.fpc - this.n)/this.fpc
      }
      
      xs<-x[this.stratum,,drop=FALSE]
      
      ## stratum with only 1 cluster leads to undefined variance
      if (this.n==1){
          this.df<-1
          lonely.psu<-match.arg(lonely.psu, c("remove","adjust","fail","certainty"))
          if (lonely.psu=="fail")
              stop("Stratum ",s, " has only one sampling unit.")
          else if (lonely.psu!="certainty")
              warning("Stratum ",s, " has only one sampling unit.")
          if (lonely.psu=="adjust")
            xs<-strata.means[match(s,ss),,drop=FALSE]
      }
      
      ## add it up
      v<-v+crossprod(xs)*this.df*this.fpc
  }
  v
}


svymean<-function(x,design, na.rm=FALSE,deff=FALSE){

  if (!inherits(design,"survey.design"))
    stop("design is not a survey design")
  
  if (inherits(x,"formula"))
    x<-model.frame(x,design$variables,na.action=na.pass)
  else if(typeof(x) %in% c("expression","symbol"))
    x<-eval(x, design$variables)
  
  x<-as.matrix(x)
  
  if (na.rm){
    nas<-rowSums(is.na(x))
            design<-design[nas==0,]
    x<-x[nas==0,,drop=FALSE]
  }
  
  pweights<-1/design$prob
  psum<-sum(pweights)
  average<-colSums(x*pweights/psum)
  x<-sweep(x,2,average)
  v<-svyCprod(x*pweights/psum,design$strata,design$cluster[[1]], design$fpc, design$nPSU)
  attr(average,"var")<-v
  attr(average,"statistic")<-"mean"
  class(average)<-"svystat"
  if (deff){
    vsrs<-svyvar(x,design,na.rm=na.rm)/NROW(design$cluster)
    attr(average, "deff")<-v/vsrs
  }
  
  return(average)
}


print.svystat<-function(x,...){
  m<-cbind(x,sqrt(diag(attr(x,"var"))))
  deff<-attr(x,"deff")
  if (is.null(deff)){
    colnames(m)<-c(attr(x,"statistic"),"SE")
  } else {
    m<-cbind(m,diag(as.matrix(deff)))
    colnames(m)<-c(attr(x,"statistic"),"SE","DEff")
  }
  printCoefmat(m)
}

svytotal<-function(x,design, na.rm=FALSE, deff=FALSE){

  if (!inherits(design,"survey.design"))
    stop("design is not a survey design")
  
  if (inherits(x,"formula"))
      x<-model.frame(x,design$variables,na.action=na.pass)
  else if(typeof(x) %in% c("expression","symbol"))
      x<-eval(x, design$variables)

  x<-as.matrix(x)
  
  if (na.rm){
    nas<-rowSums(is.na(x))
    design<-design[nas==0,]
    x<-x[nas==0,,drop=FALSE]
  }

  N<-sum(1/design$prob)
  m <- svymean(x, design, na.rm=na.rm)
  total<-m*N
  attr(total, "var")<-v<-svyCprod(x/design$prob,design$strata, design$cluster[[1]], design$fpc, design$nPSU)
  attr(total,"statistic")<-"total"
  if (deff){
    vsrs<-svyvar(x,design)*sum(weights(design)^2)
    attr(total,"deff")<-v/vsrs
  }
  return(total)
}

svyvar<-function(x, design, na.rm=FALSE){
    
	if (inherits(x,"formula"))
            x<-model.frame(x,design$variables,na.action=na.pass)
	else if(typeof(x) %in% c("expression","symbol"))
            x<-eval(x, design$variables)
        
	xbar<-svymean(x,design, na.rm=na.rm)
	if(NCOL(x)==1) {
            x<-x-xbar
            v<-svymean(x*x,design, na.rm=na.rm)
            attr(v,"statistic")<-"variance"
            return(v)
	}
	x<-t(t(x)-xbar)
	p<-NCOL(x)
	n<-NROW(x)
	a<-matrix(rep(x,p),ncol=p*p)
	b<-x[,rep(1:p,each=p)]
	v<-svymean(a*b,design, na.rm=na.rm)
	v<-matrix(v,ncol=p)
        attr(v,"statistic")<-"variance"
    }

svyquantile<-function(x,design,quantiles,alpha=0.05,ci=FALSE, method="linear",f=1){
    if (inherits(x,"formula"))
		x<-model.frame(x,design$variables)
    else if(typeof(x) %in% c("expression","symbol"))
        x<-eval(x, design$variables)
    
    w<-weights(design)
    
    computeQuantiles<-function(xx){
      oo<-order(xx)
      cum.w<-cumsum(w[oo])/sum(w)
      cdf<-approxfun(cum.w,xx[oo],method="linear",f=1,
                     yleft=min(xx),yright=max(xx)) 
      cdf(quantiles)
    }
      
    computeCI<-function(xx,p){
    
      U<-function(theta){ ((xx>theta)-(1-p))}
      
      scoretest<-function(theta,qlimit){
        umean<-svymean(U(theta),design)
        umean/sqrt(attr(umean,"var"))-qlimit
      }

      iqr<-IQR(xx)
      lower<-min(xx)+iqr/100
      upper<-max(xx)-iqr/100
      tol<-1/(100*sqrt(nrow(design)))
      c(uniroot(scoretest,interval=c(lower,upper),qlimit=qnorm(alpha/2,lower.tail=FALSE),tol=tol)$root,
        uniroot(scoretest,interval=c(lower,upper),qlimit=qnorm(alpha/2,lower.tail=TRUE),tol=tol)$root)
    }

    if (!is.null(dim(x)))
      rval<-t(matrix(apply(x,2,computeQuantiles),nrow=length(quantiles),
                   dimnames=list(as.character(round(quantiles,2)),colnames(x))))
    else
      rval<-computeQuantiles(x)

    if (!ci) return(rval)

    if (!is.null(dim(x)))
      cis<-array(apply(x,2,function(xx) sapply(quantiles,function(qq) computeCI(xx,qq))),
                 dim=c(2,length(quantiles),ncol(x)),
                 dimnames=list(c("(lower","upper)"),
                   as.character(round(quantiles,2)),
                   colnames(x)))
    else
      cis<-sapply(quantiles, function(qq) computeCI(x,qq))

    
    list(quantiles=rval,CIs=cis)
  
    
  }

svyratio<-function(numerator, denominator, design){

    if (inherits(numerator,"formula"))
		numerator<-model.frame(numerator,design$variables)
    else if(typeof(numerator) %in% c("expression","symbol"))
        numerator<-eval(numerator, design$variables)
    if (inherits(denominator,"formula"))
		denominator<-model.frame(denominator,design$variables)
    else if(typeof(denominator) %in% c("expression","symbol"))
        denominator<-eval(denominator, design$variables)

    nn<-NCOL(numerator)
    nd<-NCOL(denominator)

    all<-cbind(numerator,denominator)
    allstats<-svytotal(all,design) 
    rval<-list(ratio=outer(allstats[1:nn],allstats[nn+1:nd],"/"))


    vars<-matrix(ncol=nd,nrow=nn)
    for(i in 1:nn){
      for(j in 1:nd){
        r<-(numerator[,i]-rval$ratio[i,j]*denominator[,j])/sum(denominator[,j]/design$prob)
        vars[i,j]<-svyCprod(r*1/design$prob, design$strata, design$cluster[[1]], design$fpc, design$nPSU)
      }
    }
    colnames(vars)<-names(denominator)
    rownames(vars)<-names(numerator)
    rval$var<-vars
    rval$call<-sys.call()
    class(rval)<-"svyratio"
    rval
    
  }

print.svyratio<-function(x,...){
  cat("Ratio estimator: ")
  print(x$call)
  cat("Ratios=\n")
  print(x$ratio)
  cat("SEs=\n")
  print(sqrt(x$var))
  invisible(NULL)
}

predict.svyratio<-function(object, total, se=TRUE,...){
  if (se)
    return(list(total=object$ratio*total,se=sqrt(object$var)*total))
  else
    return(object$ratio*total)
}

svytable<-function(formula, design, Ntotal=design$fpc, round=FALSE){

  if (!inherits(design,"survey.design")) stop("design must be a survey design")
    weights<-1/design$prob
  
   ## unstratified or unadjusted.
   if (is.null(Ntotal) || length(Ntotal)==1){
       if (length(formula)==3)
           tblcall<-bquote(xtabs(I(weights*.(formula[[2]]))~.(formula[[3]]), data=design$variables))
        else
           tblcall<-bquote(xtabs(weights~.(formula[[2]]), data=design$variables))
       tbl<-eval(tblcall)
       if (!is.null(Ntotal)) {
         if(length(formula)==3)
           tbl<-tbl/sum(Ntotal)
         else
           tbl<-tbl*sum(Ntotal)/sum(tbl)
       }
       if (round)
           tbl<-round(tbl)
       return(tbl)
   }
   ## adjusted and stratified
   if (length(formula)==3)
           tblcall<-bquote(xtabs(I(weights*.(formula[[2]]))~design$strata+.(formula[[3]]), data=design$variables))
   else
           tblcall<-bquote(xtabs(weights~design$strata+.(formula[[2]]), data=design$variables))
   tbl<-eval(tblcall)

   ss<-match(sort(unique(design$strata)), Ntotal[,1])
   dm<-dim(tbl)
   layer<-prod(dm[-1])
      tbl<-sweep(tbl,1,Ntotal[ss, 2]/apply(tbl,1,sum),"*")
   tbl<-apply(tbl, 2:length(dm), sum)
   if (round)
       tbl<-round(tbl)
   class(tbl)<-c("svytable","xtabs", "table")
   attr(tbl, "call")<-match.call()
   tbl
}

svycoxph<-function(formula,design,subset=NULL,...){
    subset<-substitute(subset)
    subset<-eval(subset,design$variables,parent.frame())
    if (!is.null(subset))
        design<-design[subset,]
    
    require(survival) || stop("Needs the survival package")
    data<-design$variables 
    
    g<-match.call()
    g$design<-NULL
    g$var<-NULL
    g$weights<-quote(.survey.prob.weights)
    g[[1]]<-quote(coxph)      
    
    ##need to rescale weights for stability 
    data$.survey.prob.weights<-(1/design$prob)/sum(1/design$prob)
    if (!all(all.vars(formula) %in% names(data))) 
        stop("all variables must be in design= argument")
    g<-with(data,eval(g))
    
    nas<-attr(model.frame(g), "na.action")
    if (length(nas))
        design<-design[-nas,]
    
    
    g$var<-svyCprod(resid(g,"dfbeta",weighted=TRUE), design$strata,
                    design$cluster[[1]], design$fpc,design$nPSU)
    
    g$naive.var<-NULL
    g$wald.test<-coef(g)%*%solve(g$var,coef(g))
    g$loglik<-c(NA,NA)
    g$rscore<-NULL
    g$score<-NA
    
    class(g)<-c("svycoxph",class(g))
    g$call<-match.call()
    g$survey.design<-design
    g
}

print.svycoxph<-function(x,...){
    print(x$survey.design, varnames=FALSE, design.summaries=FALSE,...)
    NextMethod()
}

summary.svycoxph<-function(object,...){
    print(object$survey.design,varnames=FALSE, design.summaries=FALSE,...)
    NextMethod()
}

survfit.svycoxph<-function(object,...){
    stop("No survfit method for survey models")
}
extractAIC.svycoxph<-function(fit,...){
    stop("No AIC for survey models")
}

anova.svycoxph<-function(object,...){
    stop("No anova method for survey models")
}

svyglm<-function(formula,design,subset=NULL,...){

      subset<-substitute(subset)
      subset<-eval(subset, design$variables, parent.frame())
      if (!is.null(subset))
        design<-design[subset,]
      
      data<-design$variables

      g<-match.call()
      g$design<-NULL
      g$var<-NULL
      g$weights<-quote(.survey.prob.weights)
      g[[1]]<-quote(glm)      

      ##need to rescale weights for stability in binomial
      data$.survey.prob.weights<-(1/design$prob)/sum(1/design$prob)
      if (!all(all.vars(formula) %in% names(data))) 
	stop("all variables must be in design= argument")
      g<-with(data, eval(g))

      nas<-attr(model.frame(g), "na.action")
      if (length(nas))
	design<-design[-nas,]

      g$cov.unscaled<-svy.varcoef(g,design)
      
      class(g)<-c("svyglm",class(g))
      g$call<-match.call()
      g$survey.design<-design 
      g
}

print.svyglm<-function(x,...){
  print(x$survey.design, varnames=FALSE, design.summaries=FALSE,...)
  NextMethod()

}

vcov.svyglm<-function(object,...)  object$cov.unscaled


svy.varcoef<-function(glm.object,design){
    Ainv<-summary(glm.object)$cov.unscaled
    estfun<-model.matrix(glm.object)*resid(glm.object,"working")*glm.object$weights
    B<-svyCprod(estfun,design$strata,design$cluster[[1]],design$fpc, design$nPSU)
    Ainv%*%B%*%Ainv
}

residuals.svyglm<-function(object,type = c("deviance", "pearson", "working", 
    "response", "partial"),...){
	type<-match.arg(type)
	if (type=="pearson"){
   	   y <- object$y
	   mu <- object$fitted.values
    	   wts <- object$prior.weights
           pwts<- 1/object$survey.design$prob
           pwts<- pwts/sum(pwts)
	   r<-(y - mu) * sqrt(wts/pwts)/(sqrt(object$family$variance(mu)))
	   if (is.null(object$na.action)) 
        	r
    	   else 
	        naresid(object$na.action, r)
	} else 
		NextMethod()

}

summary.svyglm<-function (object, correlation = FALSE, ...) 
{
    Qr <- object$qr
    est.disp <- TRUE
    df.r <- object$df.residual
    dispersion<-svyvar(na.omit(resid(object,"pearson")), object$survey.design)
    coef.p <- coef(object)
    covmat<-vcov(object)
    dimnames(covmat) <- list(names(coef.p), names(coef.p))
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    dn <- c("Estimate", "Std. Error")
    if (!est.disp) {
        pvalue <- 2 * pnorm(-abs(tvalue))
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", 
            "Pr(>|z|)"))
    }
    else if (df.r > 0) {
        pvalue <- 2 * pt(-abs(tvalue), df.r)
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", 
            "Pr(>|t|)"))
    }
    else {
        coef.table <- cbind(coef.p, Inf)
        dimnames(coef.table) <- list(names(coef.p), dn)
    }
    ans <- c(object[c("call", "terms", "family", "deviance", 
        "aic", "contrasts", "df.residual", "null.deviance", "df.null", 
        "iter")], list(deviance.resid = residuals(object, type = "deviance"), 
        aic = object$aic, coefficients = coef.table, dispersion = dispersion, 
        df = c(object$rank, df.r,NCOL(Qr$qr)), cov.unscaled = covmat, 
        cov.scaled = covmat))
    if (correlation) {
        dd <- sqrt(diag(covmat))
        ans$correlation <- covmat/outer(dd, dd)
    }
    ans$aliased<-is.na(object$coef)
    ans$survey.design<-list(call=object$survey.design$call)
    class(ans) <- c("summary.svyglm","summary.glm")
    return(ans)
}

print.summary.svyglm<-function(x,...){
  print(x$survey.design$call,varnames=FALSE,design.summaries=FALSE,...)
  NextMethod("print")
}

logLik.svyglm<-function(object,...){
   stop("svyglm not fitted by maximum likelihood.")
}

extractAIC.svyglm<-function(fit,...){
    stop("svyglm not fitted by maximum likelihood")
}


svymle<-function(loglike, gradient=NULL, design, formulas, start=NULL, control=list(maxit=1000), na.action="na.fail", ...){
  
 method<-if(is.null(gradient)) "Nelder-Mead" else "BFGS"

  if (!inherits(design,"survey.design")) 
	stop("design is not a survey.design")

  weights<-1/design$prob
  wtotal<-sum(weights)
  if (is.null(control$fnscale))
      control$fnscale<- -wtotal
  data<-design$variables

## Get the response variable
  nms<-names(formulas)
  if (nms[1]==""){
	if (inherits(formulas[[1]],"formula"))
	  y<-eval.parent(model.frame(formulas[[1]],data=data,na.action=na.pass))
	else
	  y<-eval(y,data,parent.frame())
	formulas[1]<-NULL
	if (NCOL(y)>1) stop("Y has more than one column")
    }   else {
  	## one formula must have response
	has.response<-sapply(formulas,length)==3
	if (sum(has.response)!=1) stop("Need a response variable")
	ff<-formulas[[which(has.response)]]
	ff[[3]]<-1
	y<-eval.parent(model.frame(ff,data=data,na.action=na.pass))
	formulas[[which(has.response)]]<-delete.response(terms(formulas[[which(has.response)]]))
        nms<-c("",nms)
  }

  if(length(which(nms==""))>1) stop("Formulas must have names")
  
  
  mf<-vector("list",length(formulas))
  for(i in 1:length(formulas)){
	mf[[i]]<-eval.parent(model.frame(formulas[[i]], data=data, na.action=na.pass))
  	if (NCOL(mf[[i]])==0) mf[[i]]<-NULL
	}
  mf<-as.data.frame(do.call("cbind",c(y,mf)))
  names(mf)[1]<-"(Response)"
  mf<-mf[,!duplicated(colnames(mf)),drop=FALSE]

  mf<-get(na.action)(mf)  
  nas<-attr(mf,"na.action")
  if (length(nas))
	design<-design[-nas,]

  Y<-mf[,1]
  mm<-lapply(formulas,model.matrix, data=mf)

  ## parameter names
  parnms<-lapply(mm,colnames)
  for(i in 1:length(parnms))
	parnms[[i]]<-paste(nms[i+1],parnms[[i]],sep=".")
  parnms<-unlist(parnms)

  # maps position in theta to model matrices
  np<-c(0,cumsum(sapply(mm,NCOL)))


  objectivefn<-function(theta,...){
     args<-vector("list",length(nms))
     args[[1]]<-Y
     for(i in 2:length(nms))
	args[[i]]<-mm[[i-1]]%*%theta[(np[i-1]+1):np[i]]
     names(args)<-nms
     args<-c(args, ...)
     sum(do.call("loglike",args)*weights)
  }

  if (is.null(gradient)) {
     grad<-NULL
  } else {  
     fnargs<-names(formals(loglike))[-1]
     grargs<-names(formals(gradient))[-1]
     if(!identical(fnargs,grargs)) stop("loglike and gradient have different arguments.")
     reorder<-na.omit(match(grargs,nms[-1]))
     ##FIXME: need to convert d/deta into d/dtheta using modelmatrix.
     grad<-function(theta,...){
       args<-vector("list",length(nms))
       args[[1]]<-Y
       for(i in 2:length(nms))
	  args[[i]]<-drop(mm[[i-1]]%*%theta[(np[i-1]+1):np[i]])
       names(args)<-nms
       args<-c(args,...)
       rval<-NULL
       tmp<-do.call("gradient",args)
       for(i in reorder){
	   rval<-c(rval, colSums(as.matrix(tmp[,i]*weights*mm[[i]])))
	}
       drop(rval)
     }
  }

  theta0<-numeric(np[length(np)])
  if (is.list(start))
      st<-do.call("c",start)
  else
      st<-start

  if (length(st)==length(theta0)) {
	theta0<-st
  } else {
	stop("starting values wrong length")
  }

  rval<-optim(theta0, objectivefn, grad,control=control,hessian=TRUE,method=method,...)
 
  if (rval$conv!=0) warning("optim did not converge")

  names(rval$par)<-parnms
  dimnames(rval$hessian)<-list(parnms,parnms)

  if (is.null(gradient)) {
	rval$invinf<-solve(-rval$hessian)
	rval$scores<-NULL
	rval$sandwich<-NULL
    }  else {
       theta<-rval$par
       args<-vector("list",length(nms))
       args[[1]]<-Y
       for(i in 2:length(nms))
	  args[[i]]<-drop(mm[[i-1]]%*%theta[(np[i-1]+1):np[i]])
       names(args)<-nms
       args<-c(args,...)
       deta<-do.call("gradient",args)
       rval$scores<-NULL
       for(i in reorder)
       	 rval$scores<-cbind(rval$scores,deta[,i]*weights*mm[[i]])

       rval$invinf<-solve(-rval$hessian)
       dimnames(rval$invinf)<-list(parnms,parnms)

       db<-rval$scores%*%rval$invinf

       rval$sandwich<-svyCprod(db,design$strata,design$psu, design$fpc, design$nPSU)
       dimnames(rval$sandwich)<-list(parnms,parnms)
     }
  rval$call<-match.call()
  rval$design<-design
  class(rval)<-"svymle"
  rval

}

coef.svymle<-function(object,...) object$par

vcov.svymle<-function(object,stderr=c("robust","model"),...) {
    stderr<-match.arg(stderr)
    if (stderr=="robust"){
	rval<-object$sandwich
	if (is.null(rval)) {
		p<-length(coef(object))
		rval<-matrix(NA,p,p)
	}
    } else {
        rval<-object$invinf*mean(1/object$design$prob)
    }
    rval
}


print.svymle<-function(x,...){
  cat("Survey-sampled mle: \n")
  print(x$call)
  cat("Coef:  \n")
  print(x$par)
}

summary.svymle<-function(object,stderr=c("robust","model"),...){
    cat("Survey-sampled mle: \n")
    print(object$call)
    stderr<-match.arg(stderr)
    tbl<-data.frame(Coef=coef(object),SE=sqrt(diag(vcov(object,stderr=stderr))))
    tbl$p.value<-format.pval(2*(1-pnorm(abs(tbl$Coef/tbl$SE))), digits=3,eps=0.001)
    print(tbl)
    print(object$design)
}



.First.lib<-function(...){
    if (is.null(getOption("survey.lonely.psu")))
        options(survey.lonely.psu="fail")
}


