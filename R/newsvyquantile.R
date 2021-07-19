svyquantile<-function(x,design,quantiles,...) UseMethod("svyquantile", design)

svyquantile.survey.design <- function (x, design, quantiles, alpha = 0.05,
                                       interval.type = c("mean", "beta","xlogit", "asin","score"),
                                       na.rm = FALSE, 
                                       ci=TRUE, se = ci,
                                       qrule=c("math","school","shahvaish","hf1","hf2","hf3","hf4","hf5","hf6","hf7","hf8","hf9"),
                                       df = NULL, ...) {
    
   if (inherits(x, "formula")) 
        x <- model.frame(x, model.frame(design), na.action = na.pass)
    else if (typeof(x) %in% c("expression", "symbol")) 
        x <- eval(x, model.frame(design, na.action = na.pass))
    if (na.rm) {
        nas <- rowSums(is.na(x))
        design <- design[nas == 0, ]
        if (length(nas) > length(design$prob)) 
            x <- x[nas == 0, , drop = FALSE]
        else x[nas > 0, ] <- 0
    }

    if(is.null(df))
        df<-degf(design)

    if (is.character(qrule)){
        qrulename<-paste("qrule",match.arg(qrule),sep="_")
        qrule<-get(qrulename, mode="function")
    }

    qcrit<-if(df==Inf) qnorm else function(...) qt(...,df=df)
    

    interval.type<-match.arg(interval.type)
    if (interval.type=="score"){
        ci_fun<-ffullerCI
    } else {
        ci_fun<-function(...) woodruffCI(...,method=interval.type)
    }

    w<-weights(design)
    
    rvals<-lapply(x,
                  function(xi){
                      r<-t(sapply(quantiles,
                                  function(p){
                                      qhat<-qrule(xi,w,p)
                                      if(ci){
                                          ci<-ci_fun(xi,qhat,p,design,qrule,alpha,df)
                                          names(ci)<-c(round(100*alpha/2,2),round(100-100*alpha/2,2))
                                          c(quantile=qhat,ci=ci,
                                            se=as.numeric(diff(ci)/(2*qcrit(1-alpha/2)))
                                            )
                                      } else {
                                          c(quantile=qhat)
                                      }
                                  }
                                  ))
                      if(!ci)
                          colnames(r)<-quantiles
                      else
                          rownames(r)<-quantiles
                      r
                  }
                  )
    attr(rvals,"hasci")<-ci
    class(rvals)<-"newsvyquantile"
    rvals
}

svyquantile.svyrep.design <- function (x, design, quantiles, alpha = 0.05,
                                       interval.type = c("mean", "beta","xlogit", "asin","quantile"),
                                       na.rm = FALSE, 
                                       ci=TRUE, se = ci,
                                       qrule=c("math","school","shahvaish","hf1","hf2","hf3","hf4","hf5","hf6","hf7","hf8","hf9"),
                                       df = NULL, return.replicates=FALSE,...) {
    interval.type <- match.arg(interval.type)
    
    if (design$type %in% c("JK1", "JKn") && interval.type == "quantile") 
        warning("Jackknife replicate weights may not give valid standard errors for quantiles")
    if (design$type %in% "other" && interval.type == "quantile") 
        message("Not all replicate weight designs give valid standard errors for quantiles.")

    if (inherits(x, "formula")) 
        x <- model.frame(x, design$variables, na.action = if (na.rm) 
            na.pass
        else na.fail)
    else if (typeof(x) %in% c("expression", "symbol")) 
        x <- eval(x, design$variables)
    if (na.rm) {
        nas <- rowSums(is.na(x))
        design <- design[nas == 0, ]
        if (length(nas) > length(design$prob)) 
            x <- x[nas == 0, , drop = FALSE]
        else x[nas > 0, ] <- 0
    }

    if (is.character(qrule)){
        qrulename<-paste("qrule",match.arg(qrule),sep="_")
        qrule<-get(qrulename, mode="function")
    }
    
    if (is.null(df)) df <- degf(design)
    if (df == Inf)  qcrit <- qnorm else qcrit <- function(...) qt(..., df = df)

    w<-weights(design,"sampling")

    if (interval.type=="quantile"){
        ci_fun<-function(...) repCI(...,return.replicates=return.replicates)
    } else {
        if(return.replicates) warning("return.replicates=TRUE only implemented for interval.type='quantile'")
        ci_fun<-woodruffCI
    }

    if ((interval.type=="quantile") && return.replicates){
              rvals<-lapply(x,
                      function(xi){ lapply(quantiles,
                                           function(p){
                                               qhat<-qrule(xi,w,p)
                                               ci<-ci_fun(xi,qhat,p,design,qrule,alpha,df)
                                               list(quantile=qhat, replicates=attr(ci,"replicates"))
                                           }
                                           )
                      }
                      )
              ests<-sapply(rvals, function(v) sapply(v, function(qi) qi$qhat))
              attr(ests, "scale") <- design$scale
              attr(ests, "rscales") <- design$rscales
              attr(ests, "mse") <- design$mse
              reps<-sapply(rvals, function(v) t(sapply(v, function(qi) qi$replicates)))
              rval<-list(quantile=ests,replicates=reps)
              
              attr(rval,"var")<-svrVar(reps, design$scale, design$rscales, mse=design$mse, coef=ests)
              class(rval)<-"svrepstat"
              return(rval)
    } else {
        rvals<-lapply(x,
                      function(xi){
                          r<- t(sapply(quantiles,
                                       function(p){
                                           qhat<-qrule(xi,w,p)
                                           if(ci){
                                               ci<-ci_fun(xi,qhat,p,design,qrule,alpha,df)
                                               names(ci)<-c(round(100*alpha/2,2),round(100-100*alpha/2,2))
                                               c(quantile=qhat,ci=ci,
                                                 se=diff(ci)/(2*qcrit(1-alpha/2))
                                                 )} else{
                                                      c(quantile=qhat)
                                                  }
                                           
                                       }
                                       ))
                          if(!ci)
                              colnames(r)<-quantiles
                          else
                              rownames(r)<-quantiles
                          r
                      }
                      )
        attr(rvals,"hasci")<-ci | se
        class(rvals)<-"newsvyquantile"
    }
    
    rvals
}

SE.newsvyquantile<-function(object,...) {
    if(!attr(object,"hasci")) stop("object does not have uncertainty estimates")
    do.call(c,lapply(object,function(ai) ai[,4]))    
}

vcov.newsvyquantile<-function(object,...) {
    if(!attr(object,"hasci")) stop("object does not have uncertainty estimates")
    r<-do.call(c,lapply(object,function(ai) ai[,4]))^2
    v<-matrix(NA,nrow=length(r),ncol=length(r))
    diag(v)<-r
    v
}

coef.newsvyquantile<-function(object,...){
    if(attr(object,"hasci")) 
        do.call(c,lapply(object,function(ai) ai[,1]))    
    else
        do.call(c,lapply(object,function(ai) ai[1,]))    
}




confint.newsvyquantile<-function(object,...){
    if(!attr(object,"hasci")) stop("object does not have uncertainty estimates")
    l<-do.call(c,lapply(object,function(ai) ai[,2]))
    u<-do.call(c,lapply(object,function(ai) ai[,3]))
    cbind(l,u)
}
    

woodruffCI<-function(x, qhat,p, design, qrule,alpha,df,method=c("mean","beta","xlogit","asin")){

    method<-match.arg(method)
    m<-svymean(x<=qhat, design)
    names(m)<-"pmed"

    pconfint<-switch(method,
                    mean=as.vector(confint(m, 1, level = 1-alpha, df = df)),
                    xlogit= {xform <- svycontrast(m, quote(log(`pmed`/(1 - `pmed`))));
                        expit(as.vector(confint(xform, 1, level = 1-alpha,  df = df)))},
                    beta={n.eff <- coef(m) * (1 - coef(m))/vcov(m);
                        rval <- coef(m)[1]
                        n.eff <- n.eff * (qt(alpha/2, nrow(design) - 1)/qt(alpha/2, degf(design)))^2
                        c(qbeta(alpha/2, n.eff * rval, n.eff * (1 - rval) +  1), qbeta(1 - alpha/2, n.eff * rval + 1, n.eff *  (1 - rval)))
                    },
                    asin={xform <- svycontrast(m, quote(asin(sqrt(`pmed`))))
                        sin(as.vector(confint(xform, 1, level = 1-alpha,  df = df)))^2
                    }
                    )
    lower<-if(is.nan(pconfint[1]) || (pconfint[1]<0)) NaN else qrule(x,weights(design,"sampling"), pconfint[1])
    upper<-if(is.nan(pconfint[2])|| (pconfint[2]>1)) NaN else qrule(x,weights(design,"sampling"), pconfint[2])
    rval<-c(lower, upper)
    names(rval)<-c(round(100*alpha/2,1),round(100*(1-alpha/2),1))
    rval
    
}



ffullerCI<-function(x, qhat, p, design, qrule, alpha, df){
        qcrit<-if(df==Inf) qnorm else function(...) qt(...,df=df)
        U <- function(theta) {
            ((x > theta) - (1 - p))
        }
        scoretest <- function(theta, qlimit) {
            umean <- svymean(U(theta), design)
            umean/SE(umean) - qlimit
        }
        iqr <- IQR(x)
        lowerT <- min(x) + iqr/100
        upperT <- max(x) - iqr/100
        tol <- 1/(100 * sqrt(nrow(design)))
        
        qlow<- uniroot(scoretest, interval = c(lowerT, upperT), qlimit = qcrit(alpha/2, lower.tail = FALSE), tol = tol)$root
        qup<-uniroot(scoretest, interval = c(lowerT, upperT), qlimit = qcrit(alpha/2, lower.tail = TRUE), tol = tol)$root

        w<-weights(design)

        c(qrule(x, w, mean(x<=qlow)), qrule(x, w, mean(x<=qup)))
        
    }


repCI<-function(x, qhat, p, design, qrule, alpha, df,return.replicates){
    qcrit<-if(df==Inf) qnorm else function(...) qt(...,df=df)

    wrep<-weights(design,"analysis")
    reps<-apply(wrep,2, function(wi) qrule(x,wi,p))
    v<-with(design, svrVar(reps, scale=scale, rscales=rscales, mse=mse,coef=qhat))

    ci<- qhat+ c(-1,1)*sqrt(v)*qcrit(1-alpha/2)

    if (return.replicates) attr(ci,"replicates")<-reps

    ci
}
