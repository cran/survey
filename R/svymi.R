
svydesign.imputationList<-function(ids, probs = NULL, strata = NULL, 
             variables = NULL, fpc = NULL, data, nest = FALSE, 
             check.strata = !nest,  weights = NULL, pps=FALSE,calibrate.formula=NULL,...){
    	designs <- lapply(data$imputations, function(d) svydesign(ids=ids, probs=probs,
              strata=strata,variables=variables,fpc=fpc,nest=nest,
              check.strata=check.strata, weights=weights,data=d,pps=pps,calibrate.formula=calibrate.formula,...))
    	rval <- list(designs=designs, call=sys.call(-1))
    	class(rval) <- "svyimputationList"
    	rval
    	}

svrepdesign.imputationList<-function(variables=NULL, repweights,weights,data,degf=NULL, mse=getOption("survey.replicates.mse"),...){
  ## dispatch on data=
  if (!is.null(variables) && !inherits(variables,"imputationList"))
    stop("'variables' must also be an 'imputationList' (or NULL)")

  
  if(!is.null(variables)){
    if (inherits(repweights,"imputationList")){
      designs <- mapply(function(v,d,r) svrepdesign(variables=v, repweights=r, weights=weights,data=NULL,degf=degf,mse=mse,...),
                        variables$imputations,data$imputations, repweights$imputations,SIMPLIFY=FALSE)
    } else {
      designs <- mapply(function(d,v) svrepdesign(variables=v, repweights=repweights, weights=weights,data=d,degf=degf,mse=mse,...),
                      data$imputations,variables$imputations,SIMPLIFY=FALSE)
    }
  }else{
    if (inherits(repweights,"imputationList")){
      designs <- mapply(function(d,r) svrepdesign(repweights=r, weights=weights,data=NULL,degf=degf,mse=mse,...),
                        data$imputations, repweights$imputations,SIMPLIFY=FALSE)
    } else {
      designs <- lapply(data$imputations, function(d) svrepdesign( repweights=repweights, weights=weights,data=d,degf=degf,mse=mse,...))
    }
  }
  rval <- list(designs=designs, call=sys.call(-1))
  class(rval) <- "svyimputationList"
  rval
}

as.svrepdesign.svyimputationList<-function(design,type=c("auto","JK1","JKn","BRR","bootstrap","subbootstrap","mrbbootstrap","Fay"),
                                           fay.rho=0, fpc=NULL, fpctype=NULL, separate.replicates=FALSE,...,compress=TRUE,
                                           mse=getOption("survey.replicates.mse")){

    if (inherits(design$design[[1]], "svyrep.design")){
        warning("already a replicate-weights design")
        return(design)
    }
    
    pwts<-weights(design$designs[[1]], "sampling")
    same.weights<-all(sapply(design$designs, function (d) isTRUE(all.equal(pwts, weights(d,"sampling")))))

    if (!same.weights && separate.replicates){
        repdesigns<-lapply(design$designs,
                           function(d) as.svrepdesign(d, type=type,fay.rho=fay.rho, fpc=fpc,fpctype=fpctype,
                                                      ..., compress=compress, mse=mse))
    } else {
        design1<-design$designs[[1]]
        repdesign1<-as.svrepdesign(design1, type=type,fay.rho=fay.rho, fpc=fpc,fpctype=fpctype,
                                   ..., compress=compress, mse=mse)

        repdesigns<-lapply(design$designs,
                           function(d) {
                               repdesign1$variables<-d$variables
                               repdesign1
                           })
    }
    rval<-list(designs=repdesigns, call=sys.call(-1))
    class(rval)<-"svyimputationList"
    rval
}



svydesign.DBimputationList<-function(ids, probs = NULL, strata = NULL, 
             variables = NULL, fpc = NULL, data, nest = FALSE, 
             check.strata = !nest,  weights = NULL, ...){
 
  design.vars<-c(all.vars(ids), all.vars(probs), all.vars(strata),all.vars(fpc), all.vars(weights))
  design.query<-paste("select", paste(design.vars,collapse=","), "from", data$imputations[1])
  design.data<-DBI::dbGetQuery(data$db$connection, design.query)

  rval<-list()
  rval$design<-svydesign(ids=ids, probs=probs, strata=strata, data=design.data,
                  fpc=fpc, variables=variables, nest=nest,check.strata=check.strata,
                  weights=weights)
  class(rval$design)<-c("DBIsvydesign", class(rval$design))
  
  rval$design$updates<-data$updates
  rval$db<-data$db
  rval$imputations<-data$imputations
  rval$variables<-NULL
  rval$call<-sys.call(-1)
    class(rval)<-"svyDBimputationList"
  rval
}

print.svyDBimputationList<-function(x,...){
  cat("DB-backed Multiple (",length(x$imputations),") imputations: ",sep="")
  print(x$call)
}

print.svyimputationList<-function(x,...){
  cat("Multiple (",length(x$designs),") imputations: ",sep="")
  print(x$call)
}

dim.svyimputationList<-function(x){
  c(dim(x$designs[[1]]),length(x$designs))
}

dimnames.svyimputationList<-function(x){
   c(dimnames(x$designs[[1]]),list(paste("imputation",1:length(x$designs))))
}

subset.svyimputationList<-function(x, subset,...){
    n<-nrow(x$designs[[1]])
    e<-substitute(subset)
    r<-eval(e,x$designs[[1]]$variables, parent.frame())
    r <- r & !is.na( r )
    x$designs[[1]]<-x$designs[[1]][r,]
    same<-TRUE
    for(i in 2:length(x$designs)){
      r1<-eval(e,x$designs[[i]]$variables, parent.frame())
      x$designs[[i]]<-x$designs[[i]][r1,]
      r1<-r1 & !is.na(r1)
      if (any(r!=r1)) {
        same<-FALSE
      }
    }
    if (!same) warning('subset differed between imputations')
    
    x$call<-sys.call(-1)
    x
  }

subset.svyDBimputationList<-function(x, subset,...,all=FALSE){
    n<-nrow(x$designs[[1]])
    e<-substitute(subset)
    df<-getvars(all.vars(e), x$db$connection, x$imputations[1],
                db.only=FALSE, updates=x$design$updates)
    r<-eval(e,df, parent.frame())
    same<-TRUE
    for(i in 2:length(x$imputations)){
      df<-getvars(all.vars(e), x$db$connection, x$imputations[i],
                  db.only=FALSE, updates=x$design$updates)
      
      r1<-eval(e,df, parent.frame())
      r1<-r1 & !is.na(r1)
      if (any(r!=r1)) {
        same<-FALSE
        if (all) r <- r & r1 else r<- r | r1
      }
    }
    if (!same) warning('subset differed between imputations')
    x$design<-x$design[r,]
    x$call<-sys.call(-1)
    x
  }

with.svyimputationList<-function (data, expr, fun, ..., multicore=getOption("survey.multicore")) {
    pf <- parent.frame()
    if (multicore && !requireNamespace("parallel",quietly=TRUE))
      multicore<-FALSE

    if (!is.null(match.call()$expr)) {
      expr <- substitute(expr)
      expr$design<-as.name(".design")
      if (multicore){
        results <- parallel::mclapply(data$designs,
                            function(.design) {
                            eval(expr, list(.design=.design),enclos=pf)
                          }
                            )
      } else{
        results <- lapply(data$designs,
                          function(.design) {
                            eval(expr, list(.design=.design),enclos=pf)
                          }
                          )
        
      }
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

withReplicates.svyimputationList<-function(design, theta,...,return.replicates=FALSE){
	if (!inherits(design$designs[[1]],"svyrep.design")) 
		stop("not a replicate-weights design")
		
   m<-match.call()
   m[[1]]<-as.name("withReplicates")	
   rval<-eval(bquote(with(design, .(m))))
   attr(rval,"call")<-m
   rval
}


with.svyDBimputationList<-function (data, expr,  ..., multicore=getOption("survey.multicore")) {
    pf <- parent.frame()
    if (!is.null(match.call()$expr)) {
      expr <- substitute(expr)
      expr$design<-as.name(".design")
      if (multicore && !requireNamespace("parallel")) multicore <-FALSE
      if (multicore){
        results<-parallel::mclapply(data$imputations,
                          function(tablename) {
                            close(data)
                            .design<-data$design
                            db<-data$db
                            db$tablename<-tablename
                            .design$db<-db
                            .design<-open(.design)
                            rval<-eval(expr, list(.design=.design),enclos=pf)
                            close(.design)
                            rval
                          }
                          )
      } else {
      results <- lapply(data$imputations,
                        function(tablename) {
                          .design<-data$design
                          db<-data$db
                          db$tablename<-tablename
                          .design$db<-db
                          eval(expr, list(.design=.design),enclos=pf)
                        }
                        )
    }
    }
    attr(results, "call") <- sys.call(-1)
    results
  }


update.svyDBimputationList<-function(object, ...){
  dots <- substitute(list(...))[-1]
  newnames <- names(dots)

  updates<-lapply(dots, function(dot){
    list(inputs=all.vars(dot),expression=dot)
  })

  if (is.null(object$design$updates))
    object$design$updates<-list(updates)
  else
    object$design$updates<-c(object$design$updates, list(updates))
  object
}

update.svyimputationList<-function(object, ...){
  dots <- substitute(list(...))[-1]
  newnames <- names(dots)
  for (i in seq(along = object$designs)) {
    for (j in seq(along = dots)) {
      object$designs[[i]]$variables[, newnames[j]] <- eval(dots[[j]], 
                           object$designs[[i]]$variables, parent.frame())
    }
  }
  object
}

close.svyDBimputationList<-function(con,...){
  dbcon<-con$db$connection
    DBI::dbDisconnect(dbcon)
  invisible(con)
}

open.svyDBimputationList<-function(con,...){
  
    dbdriver<-DBI::dbDriver(con$db$dbtype)
    con$db$connection<-DBI::dbConnect(dbdriver,dbname=con$db$dbname,...)
  
  con
}
