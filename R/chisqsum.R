
pchisqsum<-function(x,df,a,lower.tail=TRUE,
                    method=c("satterthwaite","integration")){
  
  satterthwaite<-function(a,df){
    if(any(df>1)){
      a<-rep(a,df)
    }
    tr<-mean(a)
    tr2<-mean(a^2)/(tr^2)
    
    list(scale=tr*tr2, df=length(a)/tr2)
  }
  
  chisqphi<-function(t, df=1,a=1){
      (1-(0+2i)*a*t)^(-df/2)
    }
  
  make.integrand<-function(x,DF,A){ 
    m<-length(DF)
    
    function(t){
      n<-length(t)
      tmp<-matrix(chisqphi(rep(t,each=m),rep(DF,n),rep(A,n) ),ncol=n)
      phi<-apply(tmp,2,prod)
      rval<-Im(phi*exp(-(0+1i)*t*x)/(2*pi*t))
      rval[t==0]<-x/(2*pi)
      rval
    }
    
  }
  
  
  method<-match.arg(method)
  sat<-satterthwaite(a,df)
  guess<-pchisq(x/sat$scale,sat$df,lower.tail=lower.tail)
  
  if (method=="satterthwaite")
    return(guess)

  abstol<-guess/1000
  abstol<-pmax(1e-9, abstol)
  reltol<-rep(1/1000,length(abstol))
 
  if (method=="integration"){
    
    for(i in seq(length=length(x))){
      ff<-make.integrand(x[i],df,a)
      rval<-integrate(ff,-Inf,Inf,subdivisions=1000,
                          abs.tol=abstol[i],rel.tol=reltol[i],stop.on.error=FALSE)
      if (inherits(rval, "try-error") || rval$message!="OK"){
        warning("integration failed for x=",x[i],", using Satterthwaite approximation")
      }else
      guess[i]<- if (lower.tail) 1/2-rval$value else 1/2+rval$value
    }
    return(guess)
  }  
}


##
## presumably we need something involving null and alternative covariance matrices
## to compare two models.
##
lrtdistn<-function(dev, modelcov,varu=NULL,robu=NULL,
              method=c("satterthwaite","integration")){

   if(!xor(is.null(varu), is.null(robu)))
     stop("Must specify exactly one of varu and robu")
   
   if (is.null(robu)){
     amhalf<-chol(modelcov)
     m<-amhalf%*%varu%*%t(amhalf)
   } else {
     ahalf<-chol(solve(modelcov))
     m<-ahalf%*%robu%*%t(ahalf)
   }
   e<-eigen(m,only.values=TRUE)$values
   method<-match.arg(method)
   pchisqsum(dev, df=rep(1,NCOL(m)), a=e,
             method=method, lower.tail=FALSE)
}
