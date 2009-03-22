
svydesign.character<-function (ids, probs = NULL, strata = NULL, variables = NULL, 
                               fpc = NULL, data, nest = FALSE, check.strata = !nest,
                               weights = NULL,pps=FALSE,
                               dbtype="SQLite", dbname,
                               ...) 
{

  if (dbtype == "ODBC"){
    library(RODBC)
    if (dbname=="")
      dbconn<-odbcDriverConnect(dbname,...)
    else
      dbconn<-odbcConnect(dbname,...)
  } else {
    db<-dbDriver(dbtype)
    dbconn<- dbConnect(db, dbname,...)
  }
  design.vars<-c(all.vars(ids), all.vars(probs), all.vars(strata),
                 all.vars(fpc), all.vars(weights))
  design.query<-paste("select", paste(design.vars,collapse=","), "from", data)
  if (dbtype=="ODBC")
    design.data<-sqlQuery(dbconn, design.query)
  else
    design.data<-dbGetQuery(dbconn, design.query)
    
  rval<-svydesign(ids=ids, probs=probs, strata=strata, data=design.data,
                  fpc=fpc, variables=variables, nest=nest,check.strata=check.strata,
                  weights=weights)
  rval$db<-list(dbname=dbname, tablename=data, connection=dbconn, dbtype=dbtype)
  rval$variables<-NULL
  rval$call<-sys.call(-1)
  if (dbtype=="ODBC")
    class(rval)<-c("ODBCsvydesign",class(rval))
  else
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
  design$variables<-getvars(x, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svymean",design)
}


svytotal.DBIsvydesign<-function(x, design,na.rm=FALSE,...){
  design$variables<-getvars(x, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svytotal",design)
}

svyquantile.DBIsvydesign<-function(x, design,quantiles,...){
  design$variables<-getvars(x, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svyquantile",design)
}


dropFactor<-function(mf, w){
  if(!any(w==0)) return(mf)
  dropped<-w==0
  for(i in 1:ncol(mf)) {
    if (is.factor(mf[[i]])){
      fi<-mf[[i]]
      if (all(dropped[fi==levels(fi)[1]])){
        tt<-table(fi[!dropped])
        l<-min(which(tt>0))
        levs<-levels(fi)
        mf[[i]]<-relevel(mf[[i]],ref=levs[l])
      }
    }
  }
  mf
}

svyglm.DBIsvydesign<-function(formula, design,...){
  design$variables<-dropFactor(getvars(formula, design$db$connection, design$db$tablename,updates=design$updates),
                               weights(design))
  NextMethod("svyglm",design)
}



svyplot.DBIsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  design$variables[weights(design)==0,]<-NA
  NextMethod("svyplot",design)
}


svycoplot.DBIsvydesign<-function(formula,design, style=c("hexbin","transparent"),
                            basecol="black",alpha=c(0,0.8),hexscale=c("relative","absolute"),...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename, updates=design$updates)
  design$variables[weights(design)==0,]<-NA
  NextMethod("svycoplot",design)
}

svyboxplot.DBIsvydesign<-function(formula,design, ...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  design$variables[weights(design)==0,]<-NA
  class(design)<-setdiff(class(design),"DBIsvydesign")
  svyboxplot(formula,design,...)
}


svycdf.DBIsvydesign<-function(formula,design, na.rm=TRUE, ...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svycdf",design)

}

svyolr.DBIsvydesign<-function(formula,design,...){
  design$variables<-dropFactor(getvars(formula, design$db$connection, design$db$tablename,updates=design$updates),
                               weights(design))
  NextMethod("svyolr",design)
}

svycoxph.DBIsvydesign<-function(formula,design,...){
  design$variables<-dropFactor(getvars(formula, design$db$connection, design$db$tablename,updates=design$updates),
                               weights(design))
  NextMethod("svycoxph",design)
}

svyvar.DBIsvydesign<-function(x,design,na.rm=FALSE,...){
  design$variables<-getvars(x, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svyvar",design)
}



svykm.DBIsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svykm",design)
}


svykappa.DBIsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svykappa",design)
}


svysmooth.DBIsvydesign<-function(formula,design,method=c("locpoly","quantreg"),bandwidth,quantile,df,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svysmooth",design)
}


svychisq.DBIsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svychisq",design)
}

svyratio.DBIsvydesign<-function(numerator, denominator, design,...){
  design$variables<-cbind(getvars(numerator,design$db$connection, design$db$tablename,updates=design$updates),
                          getvars(denominator,design$db$connection, design$db$tablename,updates=design$updates))
  NextMethod("svyratio",design)

}


svyby.DBIsvydesign<-function(formula, by, design,...){
  design$variables<-cbind(getvars(formula,design$db$connection, design$db$tablename,updates=design$updates),
                          getvars(by,design$db$connection, design$db$tablename,updates=design$updates))
  class(design)<-setdiff(class(design),"DBIsvydesign")
  svyby(formula,by,design,...)
}

svytable.DBIsvydesign<-function(formula,design,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("svytable",design)
}


calibrate.DBIsvydesign<-function(design,formula,...){
  design$variables<-getvars(formula, design$db$connection, design$db$tablename,updates=design$updates)
  NextMethod("calibrate",design)
}
postStratify.DBIsvydesign<-function(design, strata, population, partial = FALSE, ...) .NotYetImplemented()




subset.DBIsvydesign<-function (x, subset, ...) 
{
    e <- substitute(subset)
    x$variables<-getvars(make.formula(all.vars(e)), x$db$connection, x$db$tablename,updates=x$updates)
    r <- eval(e, x$variables, parent.frame())
    r <- r & !is.na(r)
    x <- x[r, ]
    x$call <- sys.call(-1)
    x
}



dim.DBIsvydesign<-function(x){
  w<-weights(x)
  nrow<-sum(w!=0)
   coln<-names(dbGetQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
   if (!is.null(x$updates)){
     update.names<-do.call(c, lapply(x$updates, names))
     ncol<-length(unique(c(coln,update.names)))
   } else ncol<-length(coln)
   c(nrow,ncol)
}

dimnames.DBIsvydesign<-function(x){
   rown<-rownames(x$cluster)[weights(x)!=0]
   coln<-names(dbGetQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
   if (!is.null(x$updates)){
     update.names<-do.call(c, lapply(x$updates, names))
     coln<-unique(c(coln,update.names))
   }
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

