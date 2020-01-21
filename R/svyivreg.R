
svyivreg<-function(formula, design, ...) UseMethod("svyivreg",design)

svyivreg.survey.design<-function(formula, design,...){

    .data<-model.frame(design)
    .data$.weights<-weights(design,"sampling")
    .weights<-NULL ## make CMD check happy
    estfun<-get("estfun",mode="function")
    model<- AER::ivreg(formula, data=.data, weights=.weights)

    U<-estfun.ivreg(model)/weights(design,"sampling")
    n<-NROW(U)
    infl<- U%*%model$cov.unscaled 
    v<-vcov(svytotal(infl,  design))
    
    model$invinf<-model$cov.unscaled
    model$cov.unscaled<-v
    model$df.residual<-degf(design)+1-length(coef(model))
    model$sigma<-model$sigma/sqrt(mean(weights(design,"sampling")))
    model$call<-sys.call(-1)
    class(model)<-c("svyivreg","ivreg")
    model
}


summary.svyivreg<-function(object, df = NULL, ...){
    class(object)<-"ivreg"
    summary(object, vcov.=NULL, df=df, diagnostics=FALSE,...)
}

vcov.svyivreg<-function(object,...) object$cov.unscaled

svyivreg.svyrep.design<-function(formula, design,return.replicates=FALSE,...){
    .pweights<-NULL ## make CMD check happy

    withReplicates(design, return.replicates=return.replicates,
                   function(.weights, .data){
                       .data$.pweights<-.weights
                       m<-AER::ivreg(formula,data= .data, weights=.pweights)
                       coef(m)
                   })
                      
    }


estfun.ivreg<-function (x, ...) 
{
    xmat <- model.matrix(x)
    if (any(alias <- is.na(coef(x)))) 
        xmat <- xmat[, !alias, drop = FALSE]
    wts <- weights(x)
    if (is.null(wts)) 
        wts <- 1
    res <- residuals(x)
    rval <- as.vector(res) * wts * xmat
    attr(rval, "assign") <- NULL
    attr(rval, "contrasts") <- NULL
    return(rval)
}
