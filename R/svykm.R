
svykm<-function(formula, design,...) UseMethod("svykm",design)

svykm.survey.design<-function(formula, design, ...){
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
    s<-svykm.fit(y,weights(design))
  } else {
    x<-mf[,-1]
    if (NCOL(x)>1)
      groups<-do.call(interaction,x)
    else
      groups<-as.factor(x)
    s<-lapply(levels(groups), function(g) svykm.fit(y[groups==g],weights(design)[groups==g]))
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

plot.svykm<-function(x,xlab="time",ylab="Proportion surviving",ylim=c(0,1),...){
  plot(x$time,x$surv,xlab=xlab,ylab=ylab, type="s",ylim=ylim,...)
  invisible(x)
}


lines.svykm<-function(x,xlab="time",type="s",...){
  lines(x$time,x$surv, type="s",...)
  invisible(x)
}

plot.svykmlist<-function(x, pars=NULL, ...){
  if (!is.null(pars)) pars<-as.data.frame(pars)

  if(is.null(pars))
    plot(x[[1]],...)
  else
    do.call(plot,c(list(x[[1]]),pars[1,,drop=FALSE],...))

  m<-length(x)
  if(m==1) return
  for(i in 2:m){
    if(is.null(pars))
      lines(x[[i]],...)
    else
      do.call(lines,c(list(x[[i]]),pars[i,,drop=FALSE],...))
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
