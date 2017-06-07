

psrsq<-function(object, method=c("Cox-Snell","Nagelkerke"),...){
	UseMethod("psrsq",object)
}

psrsq.glm<-function(object, method=c("Cox-Snell","Nagelkerke"),...){
	nullmodel<-update(object,.~1)
	method<-match.arg(method)
	ell0<-as.vector(logLik(nullmodel))
	ell1<-as.vector(logLik(object))
	n<-object$df.null+1
	
	mutualinf<-  -2*(ell1-ell0)/n
	r2cs<-1-exp(mutualinf)
	if (method == "Cox-Snell") 
		return(r2cs)
	scaling<-1-exp(2*ell0/n)
	r2cs/scaling
}

psrsq.svyglm<-function(object, method=c("Cox-Snell", "Nagelkerke"),...){
	method<-match.arg(method)
	if (!(object$family$family %in% c("binomial","quasibinomial","poisson","quasipoisson")))
		stop("Only implemented for discrete data")
	w<-weights(object$survey.design,"sampling")
	N<-sum(w)
	n<-sum(object$prior.weights)
	minus2ell0<-object$null.deviance*(N/n)
	minus2ell1<-object$deviance*(N/n)
	mutualinf<-(minus2ell1-minus2ell0)/N
	r2cs<-1-exp(mutualinf)
	if (method =="Cox-Snell") 
		return(r2cs)
	if (any(w<1)) warning("Weights appear to be scaled: rsquared may be wrong")
	scaling<-1-exp(-minus2ell0/N)
	r2cs/scaling
}