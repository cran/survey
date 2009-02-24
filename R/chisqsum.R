
pchisqsum<-function(x,df,a,lower.tail=TRUE,
                    method=c("satterthwaite","integration","saddlepoint")){
  
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
      if (guess[i]< 1e-5 || guess[i]> 1-1e-5)
        next ## don't even try.
      ff<-make.integrand(x[i],df,a)
      rval<-integrate(ff,-Inf,Inf,subdivisions=10000,
                          abs.tol=abstol[i],rel.tol=reltol[i],stop.on.error=FALSE)
      if (inherits(rval, "try-error") || rval$message!="OK"){
        warning("integration failed for x=",x[i],", using Satterthwaite approximation")
      }else
      guess[i]<- if (lower.tail) 1/2-rval$value else 1/2+rval$value
    }
    return(guess)
  } else if (method=="saddlepoint"){
    for(i in seq(length=length(x))){
      lambda<-rep(a,df)
      sad<-sapply(x,saddle,lambda=lambda)
      if (lower.tail) sad<-1-sad
      guess<-ifelse(is.na(sad),guess,sad)
    }
    return(guess)
  }
}

saddle<-function(x,lambda){
  d<-max(lambda)
  lambda<-lambda/d
  x<-x/d
  k0<-function(zeta) -sum(log(1-2*zeta*lambda))/2
  kprime0<-function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
  kpprime0<-function(zeta) 2*sum(lambda^2/(1-2*zeta*lambda)^2)
  n<-length(lambda)
  if (x>1.05*sum(lambda))
    hatzeta<-uniroot(function(zeta) kprime0(zeta)-x,
                     lower=-0.01,upper=1/2-1/(3*x),tol=1e-8)$root
  else
     return(NA)
  
  w<-sign(hatzeta)*sqrt(2*(hatzeta*x-k0(hatzeta)))
  v<-hatzeta*sqrt(kpprime0(hatzeta))
  if (abs(hatzeta)<1e-3)
    pnorm(w,lower.tail=FALSE)
  else
    pnorm(w+log(v/w)/w, lower.tail=FALSE)
}

