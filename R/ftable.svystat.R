

ftable.svystat<-function(x,rownames,...){

    m<-cbind(coef(x),SE(x))
    if (is.null(rownames))
        return(as.table(m))

    statname<-if (is.list(x)) attr(x[[1]],"statistic") else attr(x,"statistic")
    
    rowdim<-sapply(rownames,length)
    
    mm<-array(m,dim=c(rowdim,NCOL(m)),
             dimnames=c(as.list(rownames),
             list(c(statname,"SE"))))

    
    ftable(mm,row.vars=length(rowdim)+0:1)


}

ftable.svrepstat<-ftable.svystat

ftable.svyby<-function(x,...){
  info<-attr(x,"svyby")
  margins<-info$margins
  dimnames<-lapply(x[,margins],levels)
  dims<-sapply(dimnames,length)
  dims<-c(dims,variable=info$nstats)
  if (info$vars){
    dims<-c(dims,2)
    dimnames<-c(dimnames,
                list(sub("^statistic\\.(.*)$","\\1",info$variables)),
                list(c(info$statistic,"SE"))
                )
  } else if (info$nstats==1){
    dimnames<-c(dimnames,list(info$statistic))

  } else {
    dimnames<-c(dimnames,list(info$variables))
  }
  rval<-array(as.matrix(x[,-margins,drop=FALSE]),dim=dims,dimnames=dimnames)
  ftable(rval,row.vars=c(1,length(dim(rval))))
}
