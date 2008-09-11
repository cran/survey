## svydesign(id=~SDPPSU6, strat=~SDPSTRA6, weight=~WTPFQX6, data="set1",dbtype="SQLite",dbname="~/nhanes/imp.db", nest=TRUE)->dhanes

ODBCgetvars<-function(formula, dbconnection, tables,db.only=TRUE){
     if(is.null(formula)) return(NULL)
     if(!inherits(formula,"formula")) return(formula)
     vars<-all.vars(formula)
     if (db.only) {
     	  in.db<-vars
     	} else{
         query<-sub("@tab@",tables,"select * from @tab@ limit 1")
         oneline<-sqlQuery(dbconnection,query)
         in.db<-vars[vars %in% names(oneline)]
     }
     
     query<-paste("select",paste(in.db,collapse=", "),"from",tables)
     df<-sqlQuery(dbconnection, query, errors=TRUE)
     is.string<-sapply(df,is.character)
     if(any(is.string)) {
         for(i in which(is.string)) df[[i]]<-as.factor(df[[i]])
     }
     df
   }


print.ODBCsvydesign<-function(x,...){
  cat("ODBC-backed ")
  NextMethod()
}

summary.ODBCsvydesign<-function(object,...){
   class(object)<-c("summary.ODBCsvydesign",class(object))
   object
}

print.summary.ODBCsvydesign<-function(x,...){
   print.survey.design2(x,varnames=TRUE,design.summaries=TRUE,...)
   invisible(x)
}

close.ODBCsvydesign<-function(con,...){
  close(con$db$connection,...)
  invisible(con)
}

open.ODBCsvydesign<-function(con,...){
  db<-dbDriver(con$db$dbtype)
  con$db$connection<-odbcReConnect(con$db$connection,...)
  con
}

svymean.ODBCsvydesign<-function(x, design,...){
  design$variables<-ODBCgetvars(x, design$db$connection, design$db$tablename)
  NextMethod("svymean",design)
}


svytotal.ODBCsvydesign<-function(x, design,na.rm=FALSE,...){
  design$variables<-ODBCgetvars(x, design$db$connection, design$db$tablename)
  NextMethod("svytotal",design)
}

svyquantile.ODBCsvydesign<-function(x, design,quantiles,...){
  design$variables<-ODBCgetvars(x, design$db$connection, design$db$tablename)
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
        l<-min(which(tt>0))-1
        levs<-levels(fi)
        levels(mf[[i]])<-c(levs[-(1:l)],levs[1:l])
      }
    }
  }
  mf
}

svyglm.ODBCsvydesign<-function(formula, design,...){
  design$variables<-dropFactor(ODBCgetvars(formula, design$db$connection, design$db$tablename),weights(design))
  NextMethod("svyglm",design)
}



svyplot.ODBCsvydesign<-function(formula,design,...){
  design$variables<-ODBCgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svyplot",design)
}


svycoplot.ODBCsvydesign<-function(formula,design, style=c("hexbin","transparent"),
                            basecol="black",alpha=c(0,0.8),hexscale=c("relative","absolute"),...){
  design$variables<-ODBCgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svycoplot",design)
}

svyboxplot.ODBCsvydesign<-function(formula,design, ...){
  design$variables<-ODBCgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svyboxplot",design)

}

svycdf.ODBCsvydesign<-function(formula,design, na.rm=TRUE, ...){
  design$variables<-ODBCgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svycdf",design)

}

svyolr.ODBCsvydesign<-function(formula,design,...){
  design$variables<-dropFactor(ODBCgetvars(formula, design$db$connection, design$db$tablename), weights(design))
  NextMethod("svyolr",design)
}

svycoxph.ODBCsvydesign<-function(formula,design,...){
  design$variables<-dropFactor(ODBCgetvars(formula, design$db$connection, design$db$tablename),weights(design))
  NextMethod("svycoxph",design)
}

svyvar.ODBCsvydesign<-function(x,design,na.rm=FALSE,...){
  design$variables<-ODBCgetvars(x, design$db$connection, design$db$tablename)
  NextMethod("svyvar",design)
}



svykm.ODBCsvydesign<-function(formula,design,...){
  design$variables<-ODBCgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svykm",design)
}


svykappa.ODBCsvydesign<-function(formula,design,...){
  design$variables<-ODBCgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svykappa",design)
}


svysmooth.ODBCsvydesign<-function(formula,design,...){
  design$variables<-ODBCgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svysmooth",design)
}


svychisq.ODBCsvydesign<-function(formula,design,...){
  design$variables<-ODBCgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svychisq",design)
}

svyratio.ODBCsvydesign<-function(numerator, denominator, design,...){
  design$variables<-cbind(ODBCgetvars(numerator,design$db$connection, design$db$tablename),
                          ODBCgetvars(denominator,design$db$connection, design$db$tablename))
  NextMethod("svyratio",design)

}


svyby.ODBCsvydesign<-function(formula, by, design,...){
  design$variables<-cbind(ODBCgetvars(formula,design$db$connection, design$db$tablename),
                          ODBCgetvars(by,design$db$connection, design$db$tablename))
  class(design)<-setdiff(class(design),"ODBCsvydesign")
  svyby(formula,by,design,...)
}

svytable.ODBCsvydesign<-function(formula,design,...){
  design$variables<-ODBCgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("svytable",design)
}

calibrate.ODBCsvydesign<-function(design,formula,...){
  design$variables<-ODBCgetvars(formula, design$db$connection, design$db$tablename)
  NextMethod("calibrate",design)
}
postStratify.ODBCsvydesign<-function(design, strata, population, partial = FALSE, ...) .NotYetImplemented()

subset.ODBCsvydesign<-function (x, subset, ...) 
{
    e <- substitute(subset)
    x$variables<-ODBCgetvars(make.formula(all.vars(e)), x$db$connection, x$db$tablename)
    r <- eval(e, x$variables, parent.frame())
    r <- r & !is.na(r)
    x <- x[r, ]
    x$call <- sys.call(-1)
    x
}

update.ODBCsvydesign<-function(object,...) stop("Not yet implemented")

dim.ODBCsvydesign<-function(x,...){
   nrow<-NROW(x$prob)
   ncol<-NCOL(sqlQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
   c(nrow,ncol)
}

dimnames.ODBCsvydesign<-function(x,...){
   rown<-rownames(x$cluster)
   coln<-names(sqlQuery(x$db$conn, paste("select * from", x$db$tablename, "limit 1")))
   list(rown,coln)
}


"[.ODBCsvydesign"<-function (x, i, ..., drop = TRUE) 
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

