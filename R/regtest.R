##deviance methods not exported, used by method="LRT"
deviance.svycoxph<-function(object,...) if(length(object$ll)==2) 2 * (object$ll[1] - object$ll[2]) else 0
deviance.coxph<-function(object,...) if(length(object$loglik)==2) 2 * (object$loglik[1] - object$loglik[2]) else 0


explicit1<-function(formula){
    if (length(formula)==1) 
        return(formula==1)
    
    if (!(formula[[1]]=="+" || formula[[1]]=="*" || formula[[1]]=="/" || formula[[1]]=="^"|| formula[[1]]=="~"))
        return(FALSE)
    
    if (length(formula)==3){
        (formula[[2]]==1) || explicit1(formula[[2]]) || explicit1(formula[[3]])
    } else {
        (formula[[2]]==1) || explicit1(formula[[2]])
    }
    
}

regTermTest<-function(model, test.terms, null=NULL, df=NULL, method=c("Wald","WorkingWald","LRT"), lrt.approximation="saddlepoint"){

  method<-match.arg(method)
  
  canonicalOrder<-function(term){
    tt<-strsplit(term,":")
    tt<-lapply(tt,sort)
    sapply(tt,paste,collapse=":")
  }
  
  if(inherits(test.terms,"formula")){
      test_intercept<-explicit1(test.terms)
      test.terms<-attr(terms(test.terms),"term.labels")
  } else test_intercept<-FALSE
    
  okbeta<-!is.na(coef(model,na.rm=FALSE)) ## na.rm for svyglm
  tt<-attr(terms(model),"term.labels")
  aa<-attr(model.matrix(model),"assign")[okbeta]
  if((inherits(model,"svyloglin") || inherits(model,"svyolr"))  && attr(terms(model),"intercept")){
      aa<-aa[-1]
  }
    
  index<-which(aa %in% match(canonicalOrder(test.terms),canonicalOrder(tt)))
  if (any(is.na(index)))
    stop("Terms didn't match:",canonicalOrder(test.terms),canonicalOrder(tt))

    if (test_intercept){
        if (attr(terms(model),"intercept"))
            index<-unique(c(1,index))
        else
            stop("model does not have an intercept")
    }
    
  beta<-coef(model)[index]

  if (!is.null(null))
    beta<-beta-null
  V<-vcov(model)[index,index]

  ## this should be rewritten as methods, but that's not happening any time soon.
  if (is.null(df)){
    if (inherits(model,"svyglm"))
      df<-model$df.residual
    else if (inherits(model, "svycoxph"))
      df<-model$degf.resid
    else if (inherits(model,"lm"))
      df<-model$df.residual
    else if (inherits(model,"coxph"))
      df<-model$n-length(coef(model))
    else if (inherits(model, "MIresult"))
      df<-min(model$df[index])
    else if (inherits(model,"svyloglin"))
      df<-model$df+1-length(index)
    else if (inherits(model, "svyolr"))
      df<-model$df.residual
    else
      df<-length(resid(model))-length(coef(model))
  }

  if (method %in% c("LRT","WorkingWald")){
    if (inherits(model,"svyglm"))
      V0<-model$naive.cov
    else if (inherits(model, "svycoxph"))
      V0<-model$inv.info
    else if (inherits(model,"lm"))
      V0<-vcov(model)
    else if (inherits(model,"coxph")){
      if (is.null(model$naive.var))
        V0<-model$var
      else
        V0<-model$naive.var
    } else if (inherits(model,"svyolr")) {
      V0<-solve(model$Hess)
    } else stop("method='LRT' not supported for this model")
    V0<-V0[index,index]
    
    if (test_intercept){
        test.formula<-make.formula(c(1,test.terms))[[2]]
    } else {
        test.formula<-make.formula(test.terms)[[2]]
    }
    
    if (!("formula" %in% names(model$call)))
      names(model$call)[[2]]<-"formula"

    if (method=="LRT"){
        model0<-eval(bquote(update(.(model), .~.-(.(test.formula)),design=.(model$survey.design),subset=NULL)),
                     environment(formula(model)))  ## same missing data
        if (inherits(model,"svyglm")){
            rescale<-mean(model$prior.weights)/mean(model0$prior.weights)
        } else if (inherits(model,"svycoxph")){
            rescale<-mean(model$weights)/mean(model0$weights)
        } else rescale<-1 ##FIXME -- probably want a generic drop_and_refit function
        chisq<-deviance(model0)*rescale-deviance(model)
    } else {
        chisq<-beta%*%solve(V0)%*%beta
    }
    
    misspec<-eigen(solve(V0)%*%V, only.values=TRUE)$values
    
    if (df==Inf)
      p<-pchisqsum(chisq,rep(1,length(misspec)),misspec,method=lrt.approximation,lower.tail=FALSE)
    else
      p<-pFsum(chisq,rep(1,length(misspec)),misspec,ddf=df,method=lrt.approximation,lower.tail=FALSE)
      
    rval<-list(call=sys.call(),mcall=model$call,chisq=chisq,
               df=length(index),test.terms=test.terms, 
               p=p,lambda=misspec,ddf=df)
    if (method=="LRT")
        class(rval)<-"regTermTestLRT"
    else
        class(rval)<- "regTermTestWW"
    return(rval)
  } 
  
  
  chisq<-beta%*%solve(V)%*%beta
  if (df<Inf){
    Ftest<-chisq/length(index)
    rval<-list(call=sys.call(),mcall=model$call, Ftest=Ftest,
             df=length(index),ddf=df,test.terms=test.terms,
             p=pf(Ftest,length(index),df,lower.tail=FALSE))
  } else {
    rval<-list(call=sys.call(),mcall=model$call,chisq=chisq,
               df=length(index),test.terms=test.terms,
               p=pchisq(chisq,length(index),lower.tail=FALSE))
  }
  class(rval)<-"regTermTest"
  rval
}


print.regTermTest<-function(x,...){
  cat("Wald test for ")
  cat(x$test.terms)
  cat("\n in ")
  print(x$mcall)
  if(is.null(x$Ftest))
    cat("Chisq = ",x$chisq," on ",x$df," df: p=",format.pval(x$p),"\n")
  else
    cat("F = ",x$Ftest," on ",x$df," and ",x$ddf," df: p=",format.pval(x$p),"\n")
  invisible(x)
}

print.regTermTestLRT<-function(x,...){
  if (is.null(x$ddf) || x$ddf==Inf)
    cat("Working (Rao-Scott) LRT for ")
  else
    cat("Working (Rao-Scott+F) LRT for ")
  cat(x$test.terms)
  cat("\n in ")
  print(x$mcall)
  chisq<-x$chisq/mean(x$lambda)
  cat("Working 2logLR = ",chisq, 'p=',format.pval(x$p),"\n")
  if (length(x$lambda)>1)
    cat("(scale factors: ",signif(x$lambda/mean(x$lambda),2),")")
  else cat("df=1")
  if (!is.null(x$ddf) && is.finite(x$ddf))
    cat(";  denominator df=",x$ddf)
  cat("\n")
  invisible(x)
}

print.regTermTestWW<-function(x,...){
  if (is.null(x$ddf) || x$ddf==Inf)
    cat("Working (Rao-Scott) Wald test for ")
  else
    cat("Working (Rao-Scott+F)  for ")
  cat(x$test.terms)
  cat("\n in ")
  print(x$mcall)
  chisq<-x$chisq/mean(x$lambda)
  cat("Working Wald statistic = ",chisq, 'p=',format.pval(x$p),"\n")
  if (length(x$lambda)>1)
    cat("(scale factors: ",signif(x$lambda/mean(x$lambda),2),")")
  else cat("df=1")
  if (!is.null(x$ddf) && is.finite(x$ddf))
    cat(";  denominator df=",x$ddf)
  cat("\n")
  invisible(x)
}

svycontrast<-function(stat, contrasts,add=FALSE,...) UseMethod("svycontrast")

match.names <- function(nms,contrasts){
  l<-length(nms)
  ll<-sapply(contrasts,length)

  if (l==0) stop("No names to match")
  if (length(unlist(sapply(contrasts,names)))==0)
      return(contrasts)
  if( !all( unlist(sapply(contrasts,names)) %in% nms))
    stop("names not matched")
  
  lapply(contrasts,
         function(con) {
           r<-numeric(l)
           names(r)<-nms
           r[names(con)]<-con
           r
         })
  
}

contrast<-function(coef,var,contrasts, influence=NULL){
  nas<-is.na(var[,1])
  drop<-nas & apply(contrasts,2,function(v) all(v==0))
  if(any(drop)){
    contrasts<-contrasts[,!drop,drop=FALSE]
    coef<-coef[!drop]
    var<-var[!drop,!drop,drop=FALSE]
  }
  if (any(is.na(coef))){
    badin<-is.na(coef)
    bad<-((contrasts!=0)%*%is.na(coef))>0
    rval<-rep(NA,NROW(contrasts))
    rval[!bad]<-contrasts[!bad,!badin,drop=FALSE]%*%coef[!badin]
    v<-matrix(NA,length(rval),length(rval))
    v[!bad,!bad]<-contrasts[!bad,!badin,drop=FALSE]%*%var[!badin,!badin,drop=FALSE]%*%t(contrasts[!bad,!badin,drop=FALSE])
    dimnames(v)<-list(names(rval),names(rval))
    rval<-drop(rval)
    if(!is.null(influence)){
        attr(rval,"influence")<- influence[,!badin,drop=FALSE]%*%t(contrasts[!bad,!badin,drop=FALSE])
    }
    attr(rval, "var")<-v
  } else{
    rval<-drop(contrasts%*%coef)
    v<-contrasts%*%var%*%t(contrasts)
    dimnames(v)<-list(names(rval),names(rval))
    if(!is.null(influence)){
        attr(rval,"influence")<- influence%*%t(contrasts)
    }
    attr(rval,"var")<-v
  }
  rval
}

addQuote<-function(contrasts, original){
    ll<-as.list(original)
    names(ll)<-original
    ll<-lapply(ll, as.name)
    c(contrasts,ll)
}

addLin<-function(contrasts, original){
    id<-diag(length(original))
    dimnames(id)<-list(original,original)
    rbind(contrasts,id)
    }

svycontrast.svystat<-function(stat, contrasts,add=FALSE,...){
    if (!is.list(contrasts))
        contrasts<-list(contrast=contrasts)
    if (is.language(contrasts[[1]])){
        if(add){
            contrasts<-addQuote(contrasts,names(coef(stat)))
        }
        rval<-nlcon(contrasts,as.list(coef(stat)), vcov(stat), attr(stat,"influence"))
        class(rval)<-"svrepstat"
        attr(rval,"statistic")<-"nlcon"
        return(rval)
    }
  contrasts<-match.names(names(coef(stat)),contrasts)
    contrasts<-do.call(rbind,contrasts)
    if (add)
          contrasts<-addLin(contrasts, names(coef(stat)))
    coef<-contrast(coef(stat),vcov(stat),contrasts, attr(stat,"influence"))
  class(coef)<-"svystat"
  attr(coef,"statistic")<-"contrast"
  coef
}

svycontrast.svyolr<-function(stat, contrasts,add=FALSE,...){
  if (!is.list(contrasts))
    contrasts<-list(contrast=contrasts)
  if (is.language(contrasts[[1]])){
         if(add){
            contrasts<-addQuote(contrasts,names(coef(stat)))
         }
         rval<-nlcon(contrasts,as.list(c(coef(stat),stat$zeta)), vcov(stat))
      class(rval)<-"svystat"
      attr(rval,"statistic")<-"nlcon"
      return(rval)
  }
  contrasts <- match.names(names(coef(stat)), contrasts)
  contrasts<-do.call(rbind,contrasts)
  if (add)
          contrasts<-addLin(contrasts, names(coef(stat)))
  coef<-contrast(as.vector(as.matrix(coef(stat))),
                 vcov(stat),contrasts)
  class(coef)<-"svystat"
  attr(coef,"statistic")<-"contrast"
  coef
}




svycontrast.svrepstat<-function(stat, contrasts,add=FALSE,...){
  if (!is.list(contrasts))
    contrasts<-list(contrast=contrasts)
  if (is.language(contrasts[[1]])){
      if(add){
            contrasts<-addQuote(contrasts,names(coef(stat)))
      }
      if (is.list(stat) && !is.null(stat$replicates)){ ##replicates
        rval<-list(nlcon=nlcon(contrasts,as.list(coef(stat)),varmat=NULL))
        reps<-as.matrix(stat$replicates)
        colnames(reps)<-names(coef(stat))
        xreps<-apply(reps,1, function(repi) nlcon(datalist=as.list(repi),
                                                            exprlist=contrasts, varmat=NULL))
        rval$replicates<-if(is.matrix(xreps)) t(xreps) else as.matrix(xreps)
        attr(rval$nlcon,"var")<-svrVar(rval$replicates, scale=attr(stat$replicates,"scale"),
                                     rscales=attr(stat$replicates,"rscales"),mse=attr(stat$replicates,"mse"),
                                     coef=rval$nlcon)
      attr(rval$nlcon,"statistic")<-"nlcon"
    } else {
      rval<-nlcon(contrasts,as.list(coef(stat)), vcov(stat))
      attr(rval,"statistic")<-"nlcon"
    }
    class(rval)<-"svrepstat"
    return(rval)
  } else {
      contrasts<-match.names(names(coef(stat)), contrasts)
      contrasts<-do.call(rbind,contrasts)
      if (add)
          contrasts<-addLin(contrasts, names(coef(stat)))
     
      coef<-contrast(coef(stat), vcov(stat), contrasts)
      attr(coef,"statistic")<-"contrast"
      if (is.list(stat) && !is.null(stat$replicates)){
          coef<-list(contrast=coef,
                     replicates=tcrossprod(stat$replicates, contrasts))
      }
      class(coef)<-"svrepstat"
      coef
  }
}



nlcon<-function(exprlist, datalist, varmat, influence=NULL){
    if (!is.list(exprlist))
        exprlist<-list(contrast=exprlist)
    dexprlist<-lapply(exprlist,
                      function(expr) deriv(expr, names(datalist))[[1]])
    values<-lapply(dexprlist,
                   function(dexpr) eval(do.call(substitute, list(dexpr,datalist))))
    if (is.null(varmat))
        return(do.call(c,values))
    jac<-do.call(rbind,lapply(values,
                              function(value) attr(value,"gradient")))
    var<-jac%*%varmat%*%t(jac)
    values<-do.call(c, values)
    dimnames(var)<-list(names(values),names(values))
    attr(values, "var")<-var
    if(!is.null(influence)){
        attr(values,"influence")<-influence%*%t(jac)
    }
    values
}



svycontrast.svyglm<-svycontrast.svystat
svycontrast.svycoxph<-svycontrast.svystat
svycontrast.svrepglm<-svycontrast.svrepstat
svycontrast.svrepcoxph<-svycontrast.svrepstat


svycontrast.svyby<-svycontrast.svystat
svycontrast.default<-svycontrast.svystat


svycontrast.svyby<-function(stat, contrasts,...){

    if(!is.null(r<-attr(stat, "replicates"))){
        s<-coef(stat)
        attr(s,"var")<-vcov(stat)
        attr(s,"statistic")<-attr(stat,"svyby")$statistic
        repstat<-list(stat=s, replicates=r)
        class(repstat)<-"svrepstat"
        svycontrast(repstat, contrasts,...)
    } else NextMethod() ## default
   
}
