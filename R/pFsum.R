pFsum<-function(x,df,a,ddf=Inf,lower.tail=TRUE,...){
	if (ddf==Inf) return(pchisqsum(x,df=df,a=a,lower.tail=lower.tail,...))
	
	integrand<-function(denom){
		  pchisqsum(x*ddf/denom,df=df,a=a,lower.tail=lower.tail,...)*dchisq(denom,df=ddf)
		}
	integrate(integrand,lower=0, upper=ddf)$value+integrate(integrand,lower=ddf, upper=Inf)$value
	
	
	}
	
