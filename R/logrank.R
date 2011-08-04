svylogrank<-function(formula, design,...){
	UseMethod("svylogrank",design)
}

print.svylogrank<-function(x,...){
	m<-t(x)
	rownames(m)=""
	printCoefmat(m,has.Pvalue=TRUE,P.values=TRUE)
	invisible(NULL)
	}
	
svylogrank.survey.design2<-function(formula, design,...){
	require(survival) || stop("requires the survival package")
	tms<-delete.response(terms(formula,specials="strata"))
	findstrat<-untangle.specials(tms,"strata")
	if(length(findstrat$terms))
	   tms<-tms[-findstrat$terms]
	mf<-model.frame(tms,model.frame(design))
	if(length(mf)>1)
	   stop("Only one grouping variable allowed")
	if(!is.factor(mf[[1]]) && length(unique(mf[[1]]))>2)
	   stop("Grouping variable with more than 2 levels must be a factor")
		
	b<-coef(svycoxph(formula,design,iter=1))
	v<-vcov(svycoxph(formula,design,iter=0))
	x2<-sum(b*solve(v,b))
	rval<-c(z=b/sqrt(diag(v)), Chisq=x2, p=pchisq(x2,length(b),lower.tail=FALSE))
	class(rval)<-"svylogrank"
	rval
	}

svylogrank.twophase<-svylogrank.survey.design2
svylogrank.twophase2<-svylogrank.survey.design2

svylogrank.DBIsvydesign<-function (formula, design, ...) 
{
    design$variables <- dropFactor(getvars(formula, design$db$connection, 
        design$db$tablename, updates = design$updates, subset = design$subset), 
        weights(design))
    NextMethod("svylogrank", design)
}

svylogrank.ODBCsvydesign<-function (formula, design, ...) 
{
    design$variables <- dropFactor(getvars(formula, design$db$connection, 
        design$db$tablename, updates = design$updates), weights(design))
    NextMethod("svylogrank", design)
}	
	
svylogrank.svyrep.design<-function(formula, design,...){
	require(survival) || stop("requires the survival package")
	tms<-delete.response(terms(formula,specials="strata"))
	findstrat<-untangle.specials(tms,"strata")
	if(length(findstrat$terms))
	   tms<-tms[-findstrat$terms]
	mf<-model.frame(tms,model.frame(design))
	if(length(mf)>1)
	   stop("Only one grouping variable allowed")
	if(!is.factor(mf[[1]]) && length(unique(mf[[1]]))>2)
	   stop("Grouping variable with more than 2 levels must be a factor")
  
	rr<-withReplicates(design, function(w,df){
		  environment(formula)<-environment()
		  coef(coxph(formula,data=df,weights=w+1e-8,iter=1))
		})
	
   b<-unclass(rr)
   attr(b,"var")<-NULL
	v<-attr(rr,"var")
	x2<-sum(b*solve(v,b))
   rval<- c(z=b/sqrt(diag(as.matrix(v))), Chisq=x2, p=pchisq(x2,length(b),lower.tail=FALSE))
   class(rval)<-"svylogrank"
	rval
	}	