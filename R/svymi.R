
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
