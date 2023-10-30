

logLik.svyglm<-function(object,...){
   warning("svyglm not fitted by maximum likelihood.")
   object$deviance/-2
}

AIC.svyglm<-function(object,...,k=2,null_has_intercept=TRUE){
	if (length(list(...))){
		do.call(rbind,lapply(list(object,...),extractAIC,k=k,null_has_intercept=null_has_intercept))
    } else {
	   extractAIC(object,k=k,null_has_intercept=null_has_intercept)
    }
}
extractAIC.svyglm<-function(fit,scale,k=2,...,null_has_intercept=TRUE){
    if(is.svylm(fit)) return(extractAIC_svylm(fit,...,null_has_intercept=null_has_intercept))
    if (length(attr(terms(fit),"factors"))){
        ftest<-delete.response(formula(fit))
        if (!null_has_intercept)
            ftest<-update(ftest,.~.+`1`)
        r<-regTermTest(fit, ftest, method="LRT")
        deltabar<-mean(r$lambda)
	} else {
	    r<-list(lambda=0)
	    deltabar<-NaN
	}
	d<-fit$deviance
	c(eff.p=sum(r$lambda), AIC=d+k*sum(r$lambda),deltabar=deltabar)
}

extractAIC.svrepglm<-extractAIC.svyglm

BIC.svyglm<-function(object,...,maximal){
	if (length(list(...))){
		do.call(rbind,lapply(list(object,...),dBIC,modelM=maximal))
    } else {
	   dBIC(object,modelM=maximal)
    }
	
}
	
dBIC<-function(modela, modelM){
	pm<-modela$rank
	pM<-modelM$rank	

	if (any(!(names(coef(modela))%in% names(coef(modelM))))){
		stop("coefficients in model but not in maximal model")
		}
	index<-!(names(coef(modelM))%in% names(coef(modela)))
	n<-1+modela$df.null	
	if(any(index)){
		wald<-coef(modelM)[index]%*%solve(vcov(modelM)[index,index],coef(modelM)[index])
		detDelta<-det(solve(modelM$naive.cov[index,index,drop=FALSE],modelM$cov.unscaled[index,index,drop=FALSE]))
		dbar<-detDelta^(1/(pM-pm))
		nstar<-n/dbar	
	}else {
		wald<-0
		detDelta<-1
		dbar<-1
		nstar=NaN
		}
	c(p=pm, BIC=wald+pm*log(n)+log(detDelta)+deviance(modelM),neff=nstar)
    }

extractAIC.svrepcoxph<-function (fit, scale, k = 2, ...) .NotYetImplemented()
extractAIC.svycoxph<-function (fit, scale, k = 2, ...) 
{
    Delta<-solve(fit$inv.info, fit$var)
    deltabar <- mean(diag(Delta))
    d <- -2*fit$ll[1]
    c(eff.p = sum(diag(Delta)), AIC = d + k * sum(diag(Delta)), deltabar = deltabar)
}
AIC.svycoxph<-function (object, ..., k = 2) 
{
    if (length(list(...))) {
        do.call(rbind, lapply(list(object, ...), extractAIC, 
            k = k))
    }
    else {
        extractAIC(object, k = k)
    }
}


##  special-case the Gaussian
## dAIC=-2\max_{\beta,\sigma^2}\log\ell +2p= n\log\mathrm{RSS}-n\log n +n +n\log 2\pi +2p\bar\delta

is.svylm<-function(it) {inherits(it,"svyglm") && isTRUE(all.equal(stats::family(it),gaussian()))}

extractAIC_svylm<-function(fit,scale,k=2,...,null_has_intercept=TRUE){
    y<-fit$y
    muhat<-fit$linear.predictors
    Nhat<-sum(w<-fit$prior.weights)
    sigma2hat<-    sum((y-muhat)^2 * w) / Nhat

    minus2ellhat<- Nhat*log(sigma2hat) +Nhat +Nhat*log(2*pi)

    V0<-fit$naive.cov
    V<-vcov(fit)
    if (null_has_intercept){
        V0<-V0[-1,-1,drop=FALSE]
        V<-V[-1,-1,drop=FALSE]
    }

    ## for mu
    Delta_mu<-solve(V0*sigma2hat,V) 

    ## Now for sigma2

    Isigma2<-Nhat/(2*sigma2hat^2)
    Usigma2<- -1/(2*sigma2hat)+ (y-muhat)^2/(2*sigma2hat^2)
    Hsigma2<- sum(w*Usigma2^2)
    Deltasigma2<-Isigma2/Hsigma2
    
    ## combine
    deltabar<-mean(c(diag(Delta_mu), Deltasigma2))
    eff.p<-sum(diag(Delta_mu))+Deltasigma2
    aic <- minus2ellhat + k*eff.p
        
    c(eff.p=eff.p, AIC=aic,deltabar=deltabar)
}


