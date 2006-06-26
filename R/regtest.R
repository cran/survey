regTermTest<-function(model,test.terms, null=NULL,df=Inf){

  canonicalOrder<-function(term){
    tt<-strsplit(term,":")
    tt<-lapply(tt,sort)
    sapply(tt,paste,collapse=":")
  }
    
  if(inherits(test.terms,"formula"))
    test.terms<-attr(terms(test.terms),"term.labels")

  tt<-attr(terms(model),"term.labels")
  aa<-attr(model.matrix(model),"assign")
  if(inherits(model,"coxph") && attr(terms(model),"intercept"))
    aa<-aa[-1]
  index<-which(aa %in% match(canonicalOrder(test.terms),canonicalOrder(tt)))
  if (any(is.na(index)))
    stop("Terms didn't match:",canonicalOrder(test.terms),canonicalOrder(tt))
  
  beta<-coef(model)[index]
  if (!is.null(NULL))
    beta<-beta-null
  V<-vcov(model)[index,index]
  
  if (is.null(df)){
    if (inherits(model,"svyglm"))
      df<-model$df.residual
    else if (inherits(model, "svycoxph"))
      df<-model$degf.residual
    else if (inherits(model,"lm"))
      df<-model$df.residual
    else if (inherits(model,"coxph"))
      df<-model$n-length(coef(model))
    else if (inherits(model, "MIresult"))
      df<-min(x$df[index])
    else
      df<-length(resid(model))-length(coef(model))
  }
  
  chisq<-beta%*%solve(V)%*%beta
  if (df<Inf){
    Ftest<-chisq/length(index)
    rval<-list(call=sys.call(),mcall=model$call, Ftest=Ftest,
             df=length(index),ddf=df,test.terms=test.terms,
             p=pf(Ftest,length(index),df,lower=FALSE))
  } else {
    rval<-list(call=sys.call(),mcall=model$call,chisq=chisq,
               df=length(index),test.terms=test.terms,
               p=pchisq(chisq,length(index),lower=FALSE))
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
  invisible(NULL)
}

svycontrast<-function(stat, contrasts,...) UseMethod("svycontrast")

contrast<-function(coef,var,contrasts){
  nas<-is.na(var[,1])
  drop<-nas & apply(contrasts,2,function(v) all(v==0))
  if(any(drop)){
    contrasts<-contrasts[!drop,,drop=FALSE]
    coef<-coef[!drop]
    var<-var[!drop,!drop,drop=FALSE]
  }
  rval<-contrasts%*%coef
  attr(rval,"var")<-contrasts%*%var%*%t(contrasts)
  rval
}

svycontrast.svystat<-function(stat, contrasts,...){
  if (!is.list(contrasts))
    contrasts<-list(contrast=contrasts)
  contrasts<-do.call(rbind,contrasts)
  coef<-contrast(coef(stat),vcov(stat),contrasts)
  class(coef)<-"svystat"
  attr(coef,"statistic")<-"contrast"
  coef
}

svycontrast.svyby<-function(stat, contrasts,...){
  if (!is.list(contrasts))
    contrasts<-list(contrast=contrasts)
  contrasts<-do.call(rbind,contrasts)
  coef<-contrast(as.vector(as.matrix(coef(stat))),
                 vcov(stat),contrasts)
  class(coef)<-"svystat"
  attr(coef,"statistic")<-"contrast"
  coef
}


svycontrast.svyglm<-svycontrast.svystat
svycontrast.svycoxph<-svycontrast.svystat

svycontrast.svrepstat<-function(stat, contrasts,...){
  if (!is.list(contrasts))
    contrasts<-list(contrast=contrasts)
  contrasts<-do.call(rbind,contrasts)
  
  coef<-contrast(coef(stat),vcov(stat),contrasts)
  if (is.list(stat)){
    coef<-list(contrast=coef,
               replicates=crossprod(stat$replicates, contrasts))
  }
  class(coef)<-"svrepstat"
  attr(coef,"statistic")<-"contrast"
  coef
}

