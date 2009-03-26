
svykm<-function(formula, design, se=FALSE, ...) UseMethod("svykm",design)

svykm.survey.design<-function(formula, design,se=FALSE, ...){
  require("survival")
  if (!inherits(formula,"formula")) stop("need a formula")
  if (length(formula)!=3) stop("need a two-sided formula")
  mf<-model.frame(formula, model.frame(design), na.action=na.pass)
  mf<-na.omit(mf)
  drop<-attr(mf,"na.action")
  if (!is.null(drop)) 
    design<-design[-drop,] 
  y<-model.response(mf)
  if (!is.Surv(y) || attr(y,"type")!="right")
    stop("response must be a right-censored Surv object")      

  if (ncol(mf)==1) {
    if (se)
      s<-km.stderr(y,design)
    else 
      s<-svykm.fit(y,weights(design))
  } else {
    x<-mf[,-1]
    if (NCOL(x)>1)
      groups<-do.call(interaction,x)
    else
      groups<-as.factor(x)

    if (se){
      formula[[3]]<-1
      s<-lapply(levels(groups), function(g) svykm(formula, subset(design,groups==g),se=TRUE))
    }else{
      s<-lapply(levels(groups), function(g) svykm.fit(y[groups==g],weights(design)[groups==g]))
    }
    names(s)<-levels(groups)
    class(s)<-"svykmlist"
  }
  attr(s,"call")<-sys.call(-1)
  return(s)
}

svykm.fit<-function(y,w){
  t<-y[,"time"]
  s<-y[,"status"]
  nn<-rowsum(cbind(s,1)*w,t)
  tt<-sort(unique(t))
  N<-c(sum(w),sum(w),sum(w)-cumsum(nn[-nrow(nn),2]))
  d<-c(0,nn[,1])
  surv<-pmax(0,cumprod(1-d/N))
  rval<-list(time=c(0,tt), surv=surv)
  class(rval)<-"svykm"
  rval
}

km.stderr<-function(survobj,design){
  time<-survobj[,'time']
  status<-survobj[,'status']
  ## Brute force and ignorance: compute Y and dN as totals, use delta-method
  y<-outer(time,time,">=")[,weights(design)!=0]
  dN<-diag(status)[,weights(design)!=0]
  oo<-order(time[weights(design)!=0], -status[weights(design)!=0])
  ntimes<-length(oo)
  
  totals<-svytotal(cbind(dN[,oo],y[,oo]), design)
  y<-coef(totals)[-(1:ntimes)]
  dNbar<-coef(totals)[1:ntimes]
  
  h<-cumsum(dNbar/y)
  
  dVn<- vcov(totals)[(1:ntimes),(1:ntimes)]/outer(y,y)
  dVy <- vcov(totals)[-(1:ntimes),-(1:ntimes)]*outer(dNbar/y^2,dNbar/y^2)
  dCVny<- -vcov(totals)[(1:ntimes),-(1:ntimes)]*outer(1/y,dNbar/y^2)
  dV<-dVn+dVy+dCVny+t(dCVny)
  
  V<-numeric(ntimes)
  V[1]<-dV[1,1]
  for(i in 2:ntimes) V[i]<-V[i-1]+sum(dV[1:(i-1),i])+sum(dV[i,1:i])
  
  rval<-list(time=time[weights(design)!=0][oo],surv=exp(-h),varlog=V)
  class(rval)<-"svykm"
  rval
}


plot.svykm<-function(x,xlab="time",ylab="Proportion surviving",ylim=c(0,1),ci=NULL,lty=1,...){
  if (is.null(ci))
    ci<-!is.null(x$varlog)
  
  plot(x$time,x$surv,xlab=xlab,ylab=ylab, type="s",ylim=ylim,lty=lty,...)

  if (ci){
    lines(x$time,exp(log(x$surv)-1.96*sqrt(x$varlog)),lty=2,type="s",...)
    lines(x$time,exp(log(x$surv)+1.96*sqrt(x$varlog)),lty=2,type="s",...)
  }
  invisible(x)
}


lines.svykm<-function(x,xlab="time",type="s",ci=FALSE,lty=1,...){
  lines(x$time,x$surv, type="s",lty=lty,...)
  if (ci){
    lines(x$time,exp(log(x$surv)-1.96*sqrt(x$varlog)),lty=2,type="s",...)
    lines(x$time,exp(log(x$surv)+1.96*sqrt(x$varlog)),lty=2,type="s",...)
  }
 invisible(x)
}

plot.svykmlist<-function(x, pars=NULL, ci=FALSE,...){
  if (!is.null(pars)) pars<-as.data.frame(pars)

  if(is.null(pars))
    plot(x[[1]],ci=ci,...)
  else
    do.call(plot,c(list(x[[1]]),pars[1,,drop=FALSE],ci=ci,...))

  m<-length(x)
  if(m==1) return
  for(i in 2:m){
    if(is.null(pars))
      lines(x[[i]],ci=ci,...)
    else
      do.call(lines,c(list(x[[i]]),pars[i,,drop=FALSE],ci=ci,...))
  }
  invisible(x)
}

print.svykm<-function(x, digits=3,...,header=TRUE){
  if (header) {cat("Weighted survival curve: ")
               print(attr(x,"call"))}
  suppressWarnings({iq1<-min(which(x$surv<=0.75))
                    iq2<-min(which(x$surv<=0.5))
                    iq3<-min(which(x$surv<=0.25))})
  if (is.finite(iq1)) q1<-x$time[iq1] else q1<-Inf
  if (is.finite(iq2)) q2<-x$time[iq2] else q2<-Inf
  if (is.finite(iq3)) q3<-x$time[iq3] else q3<-Inf
 cat("Q1 =",round(q1,digits)," median =",round(q2,digits)," Q3 =",round(q3,digits),"\n")
 invisible(x)
}

print.svykmlist<-function(x, digits=3,...){
  cat("Weighted survival curves:\n")
  print(attr(x,"call"))
  for(i in 1:length(x)){
    cat(names(x)[i],": ")
    print(x[[i]],digits=digits,header=FALSE)
  }
  invisible(x)
}

quantile.svykm<-function(x, probs=c(0.75,0.5,0.25),...){
  
  iq<-sapply(probs, function(p) suppressWarnings(min(which(x$surv<=p))))
  qq<-sapply(iq, function(i) if (is.finite(i)) x$time[i] else Inf)
  names(qq)<-probs
  qq

}


confint.svykm<-function(object, parm, level=0.95,...){
  if (is.null(object$varlog)) stop("no standard errors in object")

  parm<-as.numeric(parm)
  idx<-sapply(parm, function(t) max(which(object$time<=t)))
  z<-qnorm((1-level)/2)
  ci<-exp(log(object$surv[idx])+outer(sqrt(object$varlog[idx]),c(-z,z)))
  rownames(ci)<-parm
  colnames(ci)<-format( c((1-level)/2, 1-(1-level)/2),3)
  ci
}
