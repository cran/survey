
xdesign<-function(id=NULL, strata=NULL, weights=NULL, data, fpc=NULL, adjacency=NULL,
                  overlap=c("unbiased","positive"), allow.non.binary=FALSE){
    overlap<-match.arg(overlap)

     if(is.null(id) && (overlap=="positive"))
        stop("overlap='positive' is not available for an adjacency matrix")
    if(is.null(id) && (!is.null(fpc)))
        stop("fpc is not available for an adjacency matrix")
    if (!is.null(fpc)) stop("fpc not implemented")
    if (!is.null(strata)) stop("strata not implemented")
    
    if(is.null(adjacency)){
        adjacency<-make_adjacency(id, data, overlap=overlap)
    } else if (!allow.non.binary && !is.binary(adjacency))
        stop("adjacency matrix is not binary")
    
    xdesign_adjacency(id,strata=strata, weights=weights, data=data, adjacency=adjacency)
}

is.binary<-function(object) {

    if (is.matrix(object)){
       all(object %in% c(0,1))
    } else if (inherits(object, "lMatrix") || inherits(object,"nMatrix")){
        TRUE
    } else{
        all(object@x %in% c(0,1))
    }
}

make_one_adjacency<-function(id, data){
    id<-update(id,~.-1)
    mf<-model.frame(id,data,na.action=na.fail)
    if(!is.factor(mf[[1]]) & !is.character(mf[[1]]))
        mf[[1]]<-factor(mf[[1]])
    i<-sparse.model.matrix(id, mf)
    tcrossprod(as(i,"nMatrix"))    
}

make_adjacency<-function(ids, data, overlap){
    if (length(ids)==1) warning("only one clustering dimension?")
    
    idmats<-lapply(ids, make_one_adjacency, data=data)
    switch(overlap,
           unbiased= Reduce("|", idmats),
           positive= Reduce("+", idmats)
           )
}


xdesign_adjacency<-function(id, strata=NULL, weights=NULL, data, adjacency){
    if (!is.null(strata)) stop("'strata' not available for an adjacency matrix")
    
    if (is.list(id))
        id1<-id[[1]]
    else id1<-id
    wdesign<-svydesign(ids=id1, weights=weights, strata=strata, data=data)
    
    if(!is.matrix(adjacency) && !inherits(adjacency,"Matrix"))
        stop("'adjacency' must be a matrix or Matrix")
    if(!is.binary(adjacency))
        stop("'adjacency' values must be 0 or 1")
    
    rval<-list(id=id, design=wdesign, adjacency=adjacency, call=sys.call(-1))
    class(rval)<-c("xdesign_adjacency","xdesign")
    rval
}

## Generic to allow for replicates or two-phase in the future
svyxvar<-function(design,infl,...) UseMethod("svyxvar",design)

svyxvar.xdesign<-function(design, infl,...){

    adjacency<-design$adjacency
    infl<-as.matrix(infl)
    p<-NCOL(infl)
    
    varmat<-crossprod(infl,adjacency%*%infl)/(1-mean(adjacency))
    
    varmat<-as.matrix(varmat)
    dimnames(varmat)<-list(colnames(infl),colnames(infl))
    varmat
}



tr<-function(mat) sum(diag(mat))
tr2<-function(mat) sum(mat)

## FIXME: bias correction doesn't quite agree with svymean. Modify?

degf.xdesign<-function(design,...){
    tr(design$adjacency)^2/tr2(design$adjacency)
}

svymean.xdesign<-function(x, design, na.rm=FALSE, influence=TRUE,...){
    if (!influence) stop('xdesign objects require influence=TRUE')
    rval<-svymean(x, design$design, na.rm=na.rm, influence=TRUE,...)

    varmat<-svyxvar(design, attr(rval,"influence"))
    attr(rval,"var")<-varmat
    
    rval

}

svytotal.xdesign<-function(x, design, na.rm=FALSE, influence=TRUE,...){
    if (!influence) stop('xdesign objects require influence=TRUE')
    rval<-svytotal(x, design$design, na.rm=na.rm, influence=TRUE,...)

    varmat<-svyxvar(design, attr(rval,"influence"))
    attr(rval,"var")<-varmat
    
    rval

}

svyglm.xdesign<-function (formula, design, subset = NULL, family = stats::gaussian(), 
    start = NULL, influence=TRUE,...) 
{
    if (!influence) stop('xdesign objects require influence=TRUE')
    rval<-svyglm(formula, design$design, subset=subset, family=family, start=start,  influence=TRUE,...)
    
    varmat<-svyxvar(design, attr(rval,"influence"))
    rval$cov.unscaled<-varmat
    rval$df.residual<- rval$df.residual + degf(design)-degf(design$design)
    
    rval
}



update.xdesign<-function(object,...){
    object$design<-update(object$design,...)
    object$call <- sys.call(-1)
    object
}

transform.xdesign<-function(`_data`,...){
    object$design<-update(object$design,...)
    object$call <- sys.call(-1)
    object
}

print.xdesign<- function (x, varnames = FALSE, ...) 
{
    if(is.null(x$id)){
        cat(paste0("Sparse correlated design with ",degf(x)," df:\n"))
    } else {
        k<-length(x$id)
        cat(paste0(k,"-way crossed design:\n"))
        
    }
    
    print(x$call)
    
    if (varnames) {
        cat("Data variables:\n")
            print(names(x$variables))
    }
    invisible(x)
}

## subset, svyby
## should subsetting preserve design like for surveys? Probably not.

"[.xdesign"<-function(x,i,...){
    x$design<-x$design[i,...]
    x$call <- sys.call(-1)
    if (!missing(i))
        x$adjacency<-x$adjacency[i,i]
    x
}

subset.xdesign<-function(x, subset,...) 
{
    e <- substitute(subset)
    r <- eval(e, x$design$variables, parent.frame())
    r <- r & !is.na(r)
    x <- x[r, ]
    x$call <- sys.call(-1)
    x
}


dim.xdesign<-function(x,...) dim(x$design)
dimnames.xdesign<-function(x,...) dimnames(x$design)


svyby.xdesign<-function(formula, by, design, ..., influence=TRUE){
    if (!influence) stop("xdesign objects need influence=TRUE")
        
    rval<-svyby(formula, by, design$design,...,influence=TRUE)

    varmat<-svyxvar(design, attr(rval,"influence"))
    attr(rval,"var")<-varmat
    rval$se<-sqrt(diag(varmat))
    rval$call<-sys.call(-1)

    rval
}
