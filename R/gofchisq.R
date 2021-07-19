svygofchisq<-function(formula, p, design,...){
	p<-p/sum(p)
	means<-svytotal(formula, design,...)
	rval<-chisq.test(means,p=p)
	nm<-names(coef(means))
	ncat<-length(coef(means))
	means<-svycontrast(means, list(N=rep(1,ncat)), add=TRUE)
	pN<-split(cbind(p,0*diag(p)), paste0("p_",nm))
	names(pN)<-paste0("p_",nm)
	means<-svycontrast(means, pN,add=TRUE)	
	for(i in 1:length(nm)){
		O<-as.name(nm[i])
		E<-as.name(names(pN)[i])
		expr<-list(bquote((.(O)-.(E))/sqrt(.(E))))
		names(expr)[[1]]<-paste0("X2_",O)
		means<-svycontrast(means,expr,add=TRUE)		
	}
	result<-svycontrast(means, rep(c(1,0),c(ncat,2*ncat+1)))
	lambda<-eigen(vcov(means)[1:ncat,1:ncat])$values
	tr <- mean(lambda)
    tr2 <- mean(lambda ^2)/(tr^2)
    scale = tr * tr2
    df = ncat/tr2
    rval$parameter<-c(scale=scale,df=df)
	rval$p.value<-pchisqsum(rval$statistic,rep(1,ncat), lambda,lower.tail=FALSE)
	rval$data.name<-deparse(formula)
	rval$method<-"Design-based chi-squared test for given probabilities"
	rval$lambda<-lambda
	rval
}