
svyqqmath<-function(x, design, null=qnorm, na.rm=TRUE,xlab="Expected",ylab="Observed",...){

   if (inherits(x, "formula")) 
        x <- model.frame(x, model.frame(design), na.action = na.pass)
    else if (typeof(x) %in% c("expression", "symbol")) 
        x <- eval(x, model.frame(design, na.action = na.pass))
    if (na.rm) {
        nas <- rowSums(is.na(x))
        design <- design[nas == 0, ]
        if (length(nas) > length(design$prob)) 
            x <- x[nas == 0, , drop = FALSE]
        else x[nas > 0, ] <- 0
    }
    n<-NROW(x)
    for(variable in seq_len(NCOL(x))){
        ii<-order(x[, variable])
        obsi<-x[ii, variable]
        w<-weights(design,"sampling")[ii]
        cumw<-(cumsum(w)/sum(w))*(n/(n+1))
        expi<-null(cumw)
        plot(expi,obsi,xlab=xlab,ylab=ylab,...)
    }
    invisible(NULL)
}


svyqqplot<-function(formula, design, designx=NULL, na.rm=TRUE,qrule="hf8",xlab=NULL,ylab=NULL,...){
    if (is.null(designx)){
        if (inherits(formula, "formula")) 
            x <- model.frame(formula, model.frame(design), na.action = na.pass)
        else if (typeof(x) %in% c("expression", "symbol")) 
            x <- eval(formula, model.frame(design, na.action = na.pass))
        if (na.rm) {
            nas <- rowSums(is.na(x))
            design <- design[nas == 0, ]
        if (length(nas) > length(design$prob)) 
            x <- x[nas == 0, , drop = FALSE]
        else x[nas > 0, ] <- 0
        }
        Y<-x[,1]
        X<-x[,2]
        wx<-wy<-weights(design,"sampling")
    } else {
        xform<-formula[-2]
        yform<-make.formula(formula[[2]])
        environment(yform)<-environment(formula)
        Y<- model.frame(formula, model.frame(design), na.action = na.pass)[[1]]
        wy<-weights(design,"sampling")
        X<- model.frame(formula, model.frame(designx), na.action = na.pass)[[1]]
        wx<-weights(designx,"sampling")
    }
    n<-length(Y)
    m<-length(X)
    if(is.null(xlab)) xlab<-deparse(formula[[3]])
    if(is.null(ylab)) ylab<-deparse(formula[[2]])

    if(is.character(qrule))
        qrule<-get(paste("qrule",qrule,sep="_"), mode="function")

    if (n<m){
        X<-sapply(1:n, function(i) qrule(X,wx, i/n))
    } else if (n>m){
        Y<-sapply(1:m, function(i) qrule(Y,wy, i/m))
    }
    plot(sort(X),sort(Y),xlab=xlab,ylab=ylab,...)
    
}
