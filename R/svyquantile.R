
oldsvyquantile<-function(x,design,quantiles,...) UseMethod("oldsvyquantile", design)

oldsvyquantile.survey.design<-function(x,design,quantiles,alpha=0.05,
                                    ci=FALSE, method="linear",f=1,
                                    interval.type=c("Wald","score","betaWald"),
                                    na.rm=FALSE,se=ci, ties=c("discrete","rounded"), df=NULL,...){
    if (inherits(x,"formula"))
      x<-model.frame(x ,model.frame(design), na.action=na.pass)
    else if(typeof(x) %in% c("expression","symbol"))
      x<-eval(x, model.frame(design,na.action=na.pass))
    
    if (na.rm){
        nas<-rowSums(is.na(x))
        design<-design[nas==0,]
        if (length(nas)>length(design$prob))
          x<-x[nas==0,,drop=FALSE]
        else
          x[nas>0,]<-0
      }
   

    w<-weights(design)
    epsilon<-0.1*min(w[w>0])/sum(w)


    if (is.null(df)){
      qcrit<-function(p, lower.tail=TRUE) qt(p, df=degf(design), lower.tail=lower.tail)
    } else if(df==Inf){
      qcrit <- function(p,lower.tail=TRUE) qnorm(p,lower.tail=lower.tail)
    } else {
      qcrit <- function(p,lower.tail=TRUE) qt(p,df=df,lower.tail=lower.tail)
    }


    
    computeQuantiles<-function(xx,p=quantiles){
      if (any(is.na(x))) return(NA*p)
      oo<-order(xx)
      cum.w<-cumsum(w[oo])/sum(w)
      cdf<-approxfun(cum.w,xx[oo],method=method,f=f,
                     yleft=min(xx),yright=max(xx),ties=min) 
      cdf(p)
    }
    
    computeQuantilesRounded<-function(xx,p=quantiles){
      if (any(is.na(xx))) return(NA*p)
      ww<-rowsum(w,xx,reorder=TRUE)
      xx<-sort(unique(xx))
      cum.w <- cumsum(ww)/sum(ww)
      cdf <- approxfun(cum.w, xx, method = method, f = f, 
                       yleft = min(xx), yright = max(xx),ties=min)
      cdf(p)
    }
      
    
    
    computeScoreCI<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
   
      U<-function(theta){ ((xx>theta)-(1-p))}
        
      scoretest<-function(theta,qlimit){
        umean<-svymean(U(theta),design)
        umean/sqrt(attr(umean,"var"))-qlimit
      }
      
      iqr<-IQR(xx)
      lower<-min(xx)+iqr/100
      upper<-max(xx)-iqr/100
      tol<-1/(100*sqrt(nrow(design)))
      c(uniroot(scoretest,interval=c(lower,upper),
                qlimit=qcrit(alpha/2,lower.tail=FALSE),tol=tol)$root,
        uniroot(scoretest,interval=c(lower,upper),
                qlimit=qcrit(alpha/2,lower.tail=TRUE),tol=tol)$root)
    }
    
    computePCI<-function(se,alpha,p){
      if (interval.type=="Wald"){
        p.up<-p+qcrit(alpha/2,lower.tail=FALSE)*se
        p.low<-p+qcrit(alpha/2,lower.tail=TRUE)*se
        c(p.low,p.up)
      } else if (interval.type=="betaWald"){
        n.eff <- (p*(1-p))/(se^2)
        n.eff <- n.eff * ( qt(alpha/2, nrow(design)-1)/qt(alpha/2, degf(design)) )^2
        p.up<-qbeta(1-alpha/2, n.eff*p+1, n.eff*(1-p))
        p.low<-qbeta(alpha/2,  n.eff*p, n.eff*(1-p)+1)
        c(p.low,p.up)
      }
      
    }
    
    computeWaldCI<-function(xx,p){
      if (any(is.na(xx))) return(c(NA,NA))
      theta0<-computeQuantiles(xx,p)
      U<- ((xx>theta0)-(1-p))
      wtest<-svymean(U,design)
      p.ci<-computePCI(SE(wtest),alpha,p)
      p.low<-p.ci[1]
      p.up<-p.ci[2]
      oo<-order(xx)
      cum.w<-cumsum(w[oo])/sum(w)
      approx(cum.w,xx[oo],xout=c(p.low,p.up), method=method,f=f,
             yleft=min(xx),yright=max(xx),ties=min)$y 
      
    }
    
    computeWaldCIRounded<-function(xx,p){
      if(any(is.na(xx))) return(c(NA,NA))
        theta0<-computeQuantilesRounded(xx,p)
        U<- ((xx>theta0)-(1-p))
        ww<-rowsum(w,xx, reorder=TRUE)
        uxx <- sort(unique(xx))
        wtest<-svymean(U,design)
        p.ci<-computePCI(SE(wtest),alpha,p)
        p.low<-p.ci[1]
        p.up<-p.ci[2]
        oo<-order(xx)
      cum.w<-cumsum(ww)/sum(ww)
  
        approx(cum.w,uxx,xout=c(p.low,p.up), method=method,f=f,
               yleft=min(xx),yright=max(xx),ties=min)$y 
        
      }

    ties<-match.arg(ties)
    computeQ<-switch(ties, discrete=computeQuantiles,rounded=computeQuantilesRounded)
    
    if (!is.null(dim(x)))
        rval<-t(matrix(apply(x,2,computeQ),nrow=length(quantiles),
                       dimnames=list(as.character(round(quantiles,2)),colnames(x))))
    else
      rval<-computeQ(x)
    
    if (!ci & !se) return(rval)
    
    interval.type<-match.arg(interval.type)
    
    computeCI<-switch(paste(interval.type,ties,sep="."), score.discrete=computeScoreCI,
                            score.rounded=stop("ties=\"rounded\" not available with interval.type=\"score\""),
                            Wald.rounded=computeWaldCIRounded,
                            betaWald.rounded=computeWaldCIRounded,
                            Wald.discrete=computeWaldCI,
                            betaWald.discrete=computeWaldCI)
    
    if (!is.null(dim(x)))
      cis<-array(apply(x,2,function(xx) sapply(quantiles,function(qq) computeCI(xx,qq))),
                 dim=c(2,length(quantiles),ncol(x)),
                 dimnames=list(c("(lower","upper)"),
                   as.character(round(quantiles,2)),
                   colnames(x)))
    else
      cis<-sapply(quantiles, function(qq) computeCI(x,qq))

    if (ci)
      rval<-list(quantiles=rval,CIs=cis)
    else
      rval<-list(quantiles=rval)
    
    if (is.null(dim(x)))
        ses<-(cis[2,]-cis[1,])/(2*qcrit(alpha/2,lower.tail=FALSE))
    else
        ses<-(cis[2,,]-cis[1,,])/(2*qcrit(alpha/2,lower.tail=FALSE))
    attr(rval,"SE")<-ses
    class(rval)<-c("oldsvyquantile","svyquantile")
    rval
  }



oldsvyquantile.svyrep.design<-function(x,design,quantiles,method="linear",
                                                   interval.type=c("probability","quantile"),f=1,
                                                   return.replicates=FALSE,
                                                   ties=c("discrete","rounded"),na.rm=FALSE,
                                                   alpha=0.05,df=NULL,...){

  if (!exists(".Generic",inherits=FALSE))
    .Deprecated("svyquantile")

  ties<-match.arg(ties)
  interval<-match.arg(interval.type)
  if (design$type %in% c("JK1","JKn") && interval=="quantile")
    warning("Jackknife replicate weights may not give valid standard errors for quantiles")
  if (design$type %in% "other" && interval=="quantile")
    warning("Not all replicate weight designs give valid standard errors for quantiles.")
  if (inherits(x,"formula"))
		x<-model.frame(x,design$variables,na.action=if(na.rm) na.pass else na.fail)
    else if(typeof(x) %in% c("expression","symbol"))
        x<-eval(x, design$variables)


   if (na.rm){
         nas<-rowSums(is.na(x))
       design<-design[nas==0,]
        if (length(nas)>length(design$prob))
          x<-x[nas==0,,drop=FALSE]
        else
          x[nas>0,]<-0
      }
    
  if (NROW(x)<=1){
      rval<-matrix(rep(as.matrix(x),length(quantiles)),ncol=NCOL(x),nrow=length(quantiles),byrow=TRUE)
      dimnames(rval)<-list(paste("q",round(quantiles,2),sep=""), names(x))
      if (getOption("survey.drop.replicates") && !is.null(design$selfrep) && all(design$selfrep))
          vv<-matrix(0,ncol=NCOL(x),nrow=length(quantiles))
      else
          vv<-matrix(NA,ncol=NCOL(x),nrow=length(quantiles))
      dimnames(vv)<-list(paste("q",round(quantiles,2),sep=""), names(x))
      attr(rval,"var")<-vv
      attr(rval,"statistic")<-quantiles
      if (return.replicates)
          rval<-list(mean=rval,replicates=NULL)
      class(rval)<-"svrepstat"
      return(rval)
  }

    if (is.null(df))
        df<-degf(design)
    if (df==Inf)
        qcrit<-qnorm
    else
        qcrit <-function(...) qt(...,df=df)
    
  
    w<-weights(design,"analysis")


    if (interval=="quantile"){
      ## interval on quantile scale
      if (ties=="discrete")
        computeQuantiles<-function(xx){
          oo<-order(xx)
          
          ws<-weights(design,"sampling")
          epsilon<-0.1*min(ws[ws>0])/sum(ws)

          cum.ws<-cumsum(ws[oo])/sum(ws)
 
          rval<-approx(cum.ws,xx[oo],method=method,f=f,
                       yleft=min(xx),yright=max(xx),
                       xout=quantiles,ties=min)$y
          
          cum.w<-apply(w,2,function(wi) cumsum(wi[oo])/sum(wi))
   
          qq<-apply(cum.w, 2,
                    function(cum.wi){
                        approx(cum.wi,xx[oo],method=method,f=f,
                               yleft=min(xx),yright=max(xx),
                               xout=quantiles,ties=min)$y
                    }
                    )
          if (length(quantiles)>1)
            qq<-t(qq)
          else
            qq<-as.matrix(qq)
          ##rval<-colMeans(qq)
          
          rval<-list(quantiles=rval,
                     variances=diag(as.matrix(svrVar(qq,design$scale,design$rscales,mse=design$mse,coef=rval))))
          if (return.replicates)
            rval<-c(rval, list(replicates=qq))
          rval
        } else { ##ties="rounded"
          computeQuantiles<-function(xx){
              ws<-weights(design,"sampling")
              epsilon<-0.1*min(ws[ws>0])/sum(ws)

            wws<-rowsum(ws,xx,reorder=TRUE)
            uxx<-sort(unique(xx))
            
            cum.wws<-cumsum(wws)/sum(wws)
            rval<-approx(cum.wws,uxx,method=method,f=f,
                     yleft=min(xx),yright=max(xx),
                         xout=quantiles,ties=min)$y
  
            cum.w<-apply(rowsum(w,xx,reorder=TRUE),2,function(wi) cumsum(wi)/sum(wi))
              
            qq<-apply(cum.w, 2,
                      function(cum.wi){
                          approx(cum.wi,uxx,method=method,f=f,
                                 yleft=min(xx),yright=max(xx),
                                 xout=quantiles,ties=min)$y
                      }
                      )
            if (length(quantiles)>1)
              qq<-t(qq)
            else
              qq<-as.matrix(qq)
            ##rval<-colMeans(qq)
            
            rval<-list(quantiles=rval,
                       variances=diag(as.matrix(svrVar(qq,design$scale,design$rscales,mse=design$mse,coef=rval))))
            if (return.replicates)
              rval<-c(rval, list(replicates=qq))
            rval
          }
        }
    } else {
      ## interval on probability scale, backtransformed.
      if (ties=="discrete"){
        computeQuantiles<-function(xx){
          oo<-order(xx)
          w<-weights(design,"sampling")
          epsilon<-0.1*min(w[w>0])/sum(w)

          cum.w<- cumsum(w[oo])/sum(w)
  
          Qf<-approxfun(cum.w,xx[oo],method=method,f=f,
                        yleft=min(xx),yright=max(xx),
                        ties=min)
          
          point.estimates<-Qf(quantiles)
          if(length(quantiles)==1)
            estfun<-as.numeric(xx<point.estimates)
          else
            estfun<-0+outer(xx,point.estimates,"<")
          est<-svymean(estfun,design, return.replicates=return.replicates)
          if (return.replicates)
              q.estimates<-matrix(Qf(est$replicates),nrow=NROW(est$replicates))
          zcrit<-abs(qcrit(min(alpha,1-alpha)/2))
          ci<-matrix(Qf(c(coef(est)+zcrit*SE(est), coef(est)-zcrit*SE(est))),ncol=2)
          variances<-((ci[,1]-ci[,2])/2/zcrit)^2
          rval<-list(quantiles=point.estimates,
                     variances=variances)
          if (return.replicates)
            rval<-c(rval, list(replicates=q.estimates))
          rval
        }
      } else {
        ## ties=rounded
        computeQuantiles<-function(xx){
            w<-weights(design,"sampling")
            epsilon<-0.1*min(w[w>0])/sum(w)

          ww<-rowsum(w,xx,reorder=TRUE)
          uxx<-sort(unique(xx))
          cum.w<- cumsum(ww)/sum(ww)
  
          Qf<-approxfun(cum.w,uxx,method=method,f=f,
                        yleft=min(xx),yright=max(xx),
                        ties=min)
          
          point.estimates<-Qf(quantiles)
          if(length(quantiles)==1)
            estfun<-as.numeric(xx<point.estimates)
          else
                estfun<-0+outer(xx,point.estimates,"<")
          est<-svymean(estfun, design, return.replicates=return.replicates)
          if (return.replicates)
              q.estimates<-matrix(Qf(est$replicates),nrow=NROW(est$replicates))
          zcrit<-abs(qcrit(min(alpha,1-alpha)/2))
          ci<-matrix(Qf(c(coef(est)+zcrit*SE(est), coef(est)-zcrit*SE(est))),ncol=2)
          variances<-((ci[,1]-ci[,2])/2/zcrit)^2
          rval<-list(quantiles=point.estimates,
                     variances=variances)
          if (return.replicates)
            rval<-c(rval, list(replicates=q.estimates))
          rval
        }
        
      }
    }

  if (!is.null(dim(x)))
    results<-apply(x,2,computeQuantiles)
  else
    results<-computeQuantiles(x)
  
  rval<-matrix(sapply(results,"[[","quantiles"),ncol=NCOL(x),nrow=length(quantiles),
               dimnames=list(paste("q",round(quantiles,2),sep=""), names(x)))
  vv<-matrix(sapply(results,"[[","variances"),ncol=NCOL(x),nrow=length(quantiles),
             dimnames=list(paste("q",round(quantiles,2),sep=""), names(x)))
  attr(rval,"var")<-vv
  attr(rval, "statistic")<-"quantiles"
  if (return.replicates) {
    reps<-do.call(cbind,lapply(results,"[[","replicates"))
    attr(reps,"scale")<-design$scale
    attr(reps,"rscales")<-design$rscales
    attr(reps,"mse")<-design$mse
    rval<-list(mean=rval, replicates=reps)
  }
  class(rval)<-"svrepstat"
  rval
  
}


SE.oldsvyquantile<-function(object,...){
    attr(object,"SE")
}

vcov.oldsvyquantile<-function(object,...){
  se<-SE(object)
  if (is.null(se)) stop("no uncertainty information present")
  v<-matrix(NA,length(se),length(se))
  warning("Only diagonal of vcov() available")
  diag(v)<-se
  v
}

coef.oldsvyquantile<-function(object,...){
  rval<-as.vector(object$quantiles)
  if(ncol(object$quantiles)==1)
    names(rval)<-rownames(object$quantiles)
  else if (nrow(object$quantiles)==1)
    names(rval)<-colnames(object$quantiles)
  else names(rval)<-t(outer(colnames(object$quantiles),
                            rownames(object$quantiles),
                            paste,sep=":"))
  rval
}

print.oldsvyquantile<-function(x,...){
    print(list(quantiles=x$quantiles, CIs=x$CIs))
}



confint.oldsvyquantile<-function(object,parm=NULL,level=NULL,...){
  if (!is.null(level)) stop("need to re-run svyquantile to specify level")
  ci<-t(matrix(as.vector(object$CIs),nrow=2))
  colnames(ci)<-dimnames(object$CIs)[[1]]
  rownames(ci)<-outer(dimnames(object$CIs)[[2]],
                      dimnames(object$CIs)[[3]],paste,sep="_")
  if (is.null(parm)) 
    ci
  else 
    ci[parm,,drop=FALSE]
}
