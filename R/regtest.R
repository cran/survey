regTermTest<-function(model,test.terms, null=NULL){

        if(inherits(test.terms,"formula"))
            test.terms<-attr(terms(test.terms),"term.labels")
        
	tt<-attr(terms(model),"term.labels")
	aa<-attr(model.matrix(model),"assign")
        if(inherits(model,"coxph") && attr(terms(model),"intercept"))
          aa<-aa[-1]
	index<-which(aa %in% match(test.terms,tt))
        
	beta<-coef(model)[index]
        if (!is.null(NULL))
            beta<-beta-null
	V<-vcov(model)[index,index]

	chisq<-beta%*%solve(V)%*%beta
        rval<-list(call=sys.call(),mcall=model$call,chisq=chisq,
                   df=length(index),test.terms=test.terms,
                   p=pchisq(chisq,length(index),lower=FALSE))
        class(rval)<-"regTermTest"
        rval
}

print.regTermTest<-function(x,...){
       cat("Wald test for ")
       cat(x$test.terms)
       cat("\n in ")
       print(x$mcall)
       cat("Chisq = ",x$chisq," on ",x$df," df: p=",format.pval(x$p),"\n")
       invisible(NULL)
}
