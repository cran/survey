##
## Tests for contingency tables
##


svychisq<-function(formula, design, statistic=c("F","Chisq")){
  if (ncol(attr(terms(formula),"factors"))>2)
    stop("Only 2-d tables at the moment")
  statistic<-match.arg(statistic)
  
  rows<-formula[[2]][[2]]
  cols<-formula[[2]][[3]]
  nr<-length(unique(design$variables[,as.character(rows)]))
  
  fsat<-eval(bquote(~interaction(factor(.(rows)),factor(.(cols)))-1))
  
  mm<-model.matrix(fsat,model.frame(fsat, design$variables))

  mean2<-svymean(mm,design)
  nc<-length(mean2)/nr
  N<-nrow(mm)
  
  mf1<-expand.grid(rows=1:nr,cols=1:nc)
  X1<-model.matrix(~factor(rows)+factor(cols),mf1)
  X12<-model.matrix(~factor(rows)*factor(cols),mf1)
  Cmat<-qr.resid(qr(X1),X12[,-(1:(nr+nc-1)),drop=FALSE])
  

  Dmat <- diag(mean2)
  iDmat<- diag(ifelse(mean2==0,0,1/mean2))
  Vsrs <- (Dmat - outer(mean2,mean2))/N

  V <- attr(mean2,"var")

  denom<- t(Cmat) %*% (iDmat/N) %*% Cmat
  numr<-t(Cmat)%*% iDmat %*% V %*% iDmat %*% Cmat

  Delta<-solve(denom,numr)
  
  d0<- sum(diag(Delta))^2/(sum(diag(Delta%*%Delta)))
  
  warn<-options(warn=-1) ## turn off the small-cell count warning.
  pearson<- chisq.test(svytable(formula,design,Ntotal=N),
                       correct=FALSE)
  options(warn)
  
  if (match.arg(statistic)=="F"){
    pearson$statistic<-pearson$statistic/sum(diag(Delta))
    nu <- length(unique(design$cluster[,1]))-length(unique(design$strata))
    pearson$p.value<-pf(pearson$statistic, d0, d0*nu, lower.tail=FALSE)
    attr(pearson$statistic,"names")<-"F"
    pearson$parameter<-c(ndf=d0,ddf=d0*nu)
  }  else {
    pearson$p.value<-pchisq(pearson$statistic/mean(diag(Delta)),
                               df=NCOL(Delta),lower.tail=FALSE)
    pearson$parameter<-c(df=NCOL(Delta))
  }
  
  pearson$data.name<-deparse(sys.call())
  pearson$method<-"Pearson's X^2: Rao & Scott adjustment"
  pearson
  
}
