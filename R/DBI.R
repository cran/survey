## svydesign(id=~SDPPSU6, strat=~SDPSTRA6, weight=~WTPFQX6, data="set1",dbtype="SQLite",dbname="~/nhanes/imp.db", nest=TRUE)->dhanes

DBIgetvars<-function(formula, dbconnection, tables,db.only=TRUE){
     if(is.null(formula)) return(NULL)
     if(!inherits(formula,"formula")) return(formula)
     vars<-all.vars(formula)
     if (db.only) {
     	  in.db<-vars
     	} else{
         query<-sub("@tab@",tables,"select * from @tab@ limit 1")
         oneline<-dbGetQuery(dbconnection,query)
         in.db<-vars[vars %in% names(oneline)]
     }
     
     query<-paste("select",paste(in.db,collapse=", "),"from",tables)
     df<-dbGetQuery(dbconnection, query)
     is.string<-sapply(df,is.character)
     if(any(is.string)) {
         for(i in which(is.string)) df[[i]]<-as.factor(df[[i]])
     }
     df
   }


svydesign.character<-function (ids, probs = NULL, strata = NULL, variables = NULL, 
    fpc = NULL, data, nest = FALSE, check.strata = !nest, weights = NULL, dbtype="SQLite",dbname,
    ...) 
{

  db<-dbDriver(dbtype)
  dbconn<- dbConnect(db, dbname,...)
  
  design.vars<-c(all.vars(ids), all.vars(probs), all.vars(strata),all.vars(fpc), all.vars(weights))
  design.query<-paste("select", paste(design.vars,collapse=","), "from", data)
  design.data<-dbGetQuery(dbconn, design.query)

  rval<-svydesign(ids=ids, probs=probs, strata=strata, data=design.data,
                  fpc=fpc, variables=variables, nest=nest,check.strata=check.strata,
                  weights=weights)
  rval$db<-list(dbname=dbname, tablename=data, connection=dbconn, dbtype=dbtype)
  rval$data<-NULL
  class(rval)<-c("DBIsvydesign",class(rval))
  rval
}

print.DBIsvydesign<-function(x,...){
  cat("DB-backed ")
  NextMethod()
}

summary.DBIsvydesign<-function(object,...){
   class(object)<-c("summary.DBIsvydesign",class(object))
   object
}

print.summary.DBIsvydesign<-function(x,...){
   print.survey.design2(x,varnames=TRUE,design.summaries=TRUE,...)
   invisible(x)
}

close.DBIsvydesign<-function(con,...){
  dbDisconnect(con$db$connection,...)
  invisible(con)
}

open.DBIsvydesign<-function(con,...){
  db<-dbDriver(con$db$dbtype)
  con$db$connection<-dbConnect(db, dbname=con$db$dbname,...)
  con
}

svymean.DBIsvydesign<-function(x, design,...){
  design$variables<-DBIgetvars(x, design$db$connection, design$db$tablename)
  NextMethod("svymean",design)
}


svytotal.DBIsvydesign<-function(x, design,na.rm=FALSE,...){
  design$variables<-DBIgetvars(x, design$db$connection, design$db$tablename)
  NextMethod("svytotal",design)
}

svyquantile.DBIsvydesign<-function(x, design,quantiles,...){
  design$variables<-DBIgetvars(x, design$db$connection, design$db$tablename)
  NextMethod("svyquantile",design)
}


svyglm.DBIsvydesign<-function(formula, design,...){
  design$variables<-DBIgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svyglm",design)
}



svyplot.DBIsvydesign<-function(formula,design,...){
  design$variables<-DBIgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svyplot",design)
}


svycoplot.DBIsvydesign<-function(formula,design, style=c("hexbin","transparent"),
                            basecol="black",alpha=c(0,0.8),hexscale=c("relative","absolute"),...){
  design$variables<-DBIgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svycoplot",design)
}

svyolr.DBIsvydesign<-function(formula,design,...){
  design$variables<-DBIgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svyolr",design)
}

svycoxph.DBIsvydesign<-function(formula,design,...){
  design$variables<-DBIgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svycoxph",design)
}

svyvar.DBIsvydesign<-function(x,design,na.rm=FALSE,...){
  design$variables<-DBIgetvars(x, design$db$connection, design$db$tablename)
  NextMethod("svyvar",design)
}



svykm.DBIsvydesign<-function(formula,design,...){
  design$variables<-DBIgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svykm",design)
}


svykappa.DBIsvydesign<-function(formula,design,...){
  design$variables<-DBIgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svykappa",design)
}


svysmooth.DBIsvydesign<-function(formula,design,...){
  design$variables<-DBIgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svysmooth",design)
}


svychisq.DBIsvydesign<-function(formula,design,...){
  design$variables<-DBIgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svychisq",design)
}

svyratio.DBIsvydesign<-function(numerator, denominator, design,...){
  design$variables<-cbind(DBIgetvars(numerator,design$db$connection, design$db$tablename),
                          DBIgetvars(denominator,design$db$connection, design$db$tablename))
  NextMethod("svyratio",design)

}


svyby.DBIsvydesign<-function(formula, by, design,...){
  design$variables<-cbind(DBIgetvars(formula,design$db$connection, design$db$tablename),
                          DBIgetvars(by,design$db$connection, design$db$tablename))
  class(design)<-setdiff(class(design),"DBIsvydesign")
  svyby(formula,by,design,...)
}

svytable.DBIsvydesign<-function(formula,design,...){
  design$variables<-DBIgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svytable",design)
}


calibrate.DBIsvydesign<-function(design,formula,...){
  design$variables<-DBIgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("calibrate",design)
}
postStratify.DBIsvydesign<-function(design, strata, population, partial = FALSE, ...) .NotYetImplemented()




subset.DBIsvydesign<-function (x, subset, ...) 
{
    e <- substitute(subset)
    x$variables<-DBIgetvars(make.formula(all.vars(e)), x$db$connection, x$db$tablename)
    r <- eval(e, x$variables, parent.frame())
    r <- r & !is.na(r)
    x <- x[r, ]
    x$call <- sys.call(-1)
    x
}

update.DBIsvydesign<-function(object,...) stop("Not yet implemented")

dim.DBIsvydesign<-function(x,...){
   nrow<-NROW(x$prob)
   ncol<-NCOL(dbGetQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
   c(nrow,ncol)
}

dimnames.DBIsvydesign<-function(x,...){
   rown<-rownames(x$cluster)
   coln<-names(dbGetQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
   list(rown,coln)
}


"[.DBIsvydesign"<-function (x, i, ..., drop = TRUE) 
{
  if (!missing(i)) {
    if (is.logical(i)) 
      x$prob[!i] <- Inf
    else if (is.numeric(i) && length(i)) 
      x$prob[-i] <- Inf
    else {
      tmp <- x$prob[i, ]
      x$prob <- rep(Inf, length(x$prob))
      x$prob[i, ] <- tmp
    }
    index <- is.finite(x$prob)
    psu <- !duplicated(x$cluster[index, 1])
    tt <- table(x$strata[index, 1][psu])
    if (any(tt == 1)) {
      warning(sum(tt == 1), " strata have only one PSU in this subset.")
    }
    
  }
  else {
    if (!is.null(x$variables)) 
      x$variables <- x$variables[, ..1, drop = FALSE]
  }
  x
}

