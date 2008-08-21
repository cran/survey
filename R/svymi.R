
svydesign.imputationList<-function(ids, probs = NULL, strata = NULL, 
             variables = NULL, fpc = NULL, data, nest = FALSE, 
             check.strata = !nest,  weights = NULL, ...){
    	designs <- lapply(data$imputations, function(d) svydesign(ids=ids, probs=probs,
              strata=strata,variables=variables,fpc=fpc,nest=nest,
              check.strata=check.strata, weights=weights,data=d,...))
    	rval <- list(designs=designs, call=sys.call(-1))
    	class(rval) <- "svyimputationList"
    	rval
    	}

print.svyimputationList<-function(x,...){
  cat("Multiple (",length(x$designs),") imputations: ",sep="")
  print(x$call)
}


subset.svyimputationList<-function(x, subset,...,all=FALSE){
    n<-nrow(x$designs[[1]])
    e<-substitute(subset)
    r<-eval(e,x$designs[[1]]$variables, parent.frame())
	same<-TRUE
	for(i in 2:length(x$designs)){
		r1<-eval(e,x$designs[[i]]$variables, parent.frame())
		r1<-r1 & !is.na(r1)
		if (any(r!=r1)) {
			same<-FALSE
			if (all) r <- r & r1 else r<- r | r1
		   }
		}
	if (!same) warning('subset differed between imputations')
	for(i in 1:length(x$designs)) 
		x$designs[[i]]<-x$designs[[i]][r,]
	x$call<-sys.call(-1)
	x
	}

with.svyimputationList<-function (data, expr, fun, ...) {
    pf <- parent.frame()
    if (!is.null(match.call()$expr)) {
        expr <- substitute(expr)
        expr$design<-as.name(".design")
        results <- lapply(data$designs,
                          function(.design) {
                            eval(expr, list(.design=.design),enclos=pf)
                          }
                          )
      }
    else {
      results <- lapply(data$designs, fun, ...)
    }
    if (all(sapply(results, inherits, what = "imputationResult"))) {
      class(results) <- "imputationResultList"
      results$call <- sys.call(-1)
    }
    else {
      attr(results, "call") <- sys.call(-1)
    }
    results
  }
