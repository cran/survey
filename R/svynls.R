## svynls

svynls<-function(formula, design,...){
    UseMethod("svynls", design)
}


svynls.DBIsvydesign<-function(formula, design, ...){
    design$variables <- getvars(formula, design$db$connection, 
                                design$db$tablename, updates = design$updates, subset = design$subset)
    NextMethod("svynls", design)
}


utils::globalVariables(c(".survey.prob.weight", ".survey.repwt"))

svynls.svyrep.design<-function(formula, design, start, ..., return.replicates=FALSE){
    has_vars<- intersect(all.vars(formula),colnames(design))
    dat<-model.frame(design)[,all.vars(formula)]
    meanweight<-mean(weights(design, "sampling"))
    dat$.survey.prob.weight<-weights(design, "sampling")/meanweight

    
    first<-nls(formula,dat, weights=.survey.prob.weight,start=start, ...)
    theta<-coef(first)
    
    repwts<-weights(design,"analysis")
    thetas<-matrix(0, ncol=length(theta),nrow=ncol(repwts))
    
    for(i in ncol(repwts)){
          dat$.survey.repwt<-repwts[,i]/meanweight
          model<-nls(formula,dat, .survey.repwt,start=theta, ...)
          thetas[i, ] <- coef(model)
    }

    rval<-list()
    rval$fit<-first
    
    v<-svrVar(thetas, design$scale, design$rscales, coef=theta)
    rval$naive.cov<-summary(first)$cov.unscaled
    rval$cov<-v
    rval$coef<-theta
    rval$design<-design
    rval$meanweight<-meanweight
    rval$call<-sys.call(-1)
    if (return.replicates)
        rval$replicates<-thetas
    
    class(rval)<-"svynls"
    rval
}




svynls.survey.design2<-function(formula, design, start, ..., influence=FALSE){
    
    has_vars<- intersect(all.vars(formula),colnames(design))
    dat<-model.frame(design)[,has_vars]
    meanweight<-mean(weights(design, "sampling"))
    w<-weights(design, "sampling")/meanweight
    dat$.survey.prob.weight<-w
    
    fit<-nls(formula, dat, weights=.survey.prob.weight,start=start,
             ...)
    
    v0<-summary(fit)$cov.unscaled
    theta<-coef(fit)
    
    grads<-fit$m$gradient()/sqrt(w)  ## nls scales by sqrt(weights)
    resids<-fit$m$resid()
    infl<-resids*grads%*%v0
    v<-svyrecvar(infl*w, design$cluster, design$strata, 
                 design$fpc, postStrata = design$postStrata)
    rval<-list()
    rval$fit<-fit
    rval$coef<-theta
    rval$cov<-v
    rval$naive.cov<-vcov(fit)
    rval$call<-sys.call(-1)
    rval$design<-design
    rval$meanweight<-meanweight
    if(influence)
        attr(rval,"inference")<-infl
    class(rval)<-"svynls"
    rval
}

.p.nls.convInfo<-function (x, digits, show. = getOption("show.nls.convergence", 
    TRUE)) 
{
    if (!is.null(x$convInfo)) 
        with(x$convInfo, {
            if (identical(x$call$algorithm, "port")) 
                cat("\nAlgorithm \"port\", convergence message: ", 
                  stopMessage, "\n", sep = "")
            else {
                if (!isConv || show.) {
                  cat("\nNumber of iterations", if (isConv) 
                    "to convergence:"
                  else "till stop:", finIter, "\nAchieved convergence tolerance:", 
                    format(finTol, digits = digits))
                  cat("\n")
                }
                if (!isConv) {
                  cat("Reason stopped:", stopMessage)
                  cat("\n")
                }
            }
        })
    invisible()
}


print.svynls<-function (x, digits = max(3L, getOption("digits") - 3L), ...) 
{
    cat("Nonlinear survey regression model\n")
    cat("  model: ", deparse(formula(x$fit)), "\n", sep = "")
    cat(" design: ")
    print(x$design)
    print(x$fit$m$getAllPars(), digits = digits, ...)
    cat(" ",  "weighted ", "residual sum-of-squares: ",
        format(x$fit$m$deviance()*x$meanweight, 
        digits = digits), "\n", sep = "")
    .p.nls.convInfo(x, digits = digits)
    invisible(x)
}

coef.svynls<-function(object,...) object$coef
vcov.svynls<-function(object,...) object$cov

summary.svynls <- function (object, correlation = FALSE, ...) 
{
    r <- as.vector(object$fit$m$resid())
    w <- object$fit$weights
    n <- if (!is.null(w)) 
        sum(w > 0)
    else length(r)
    param <- coef(object)
    pnames <- names(param)
    p <- length(param)
    rdf <- n - p
    resvar <- if (rdf <= 0) 
        NaN
    else deviance(object$fit)/rdf
    se <-SE(object)
    tval <- param/se
    param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(param) <- list(pnames, c("Estimate", "Std. Error", 
        "t value", "Pr(>|t|)"))
    ans <- list(formula = formula(object$fit), residuals = r, sigma = sqrt(resvar), 
                df = c(p, rdf), cov.unscaled = vcov(object)/sqrt(resvar), cov.scaled=vcov(object),
                call = object$call, 
        convInfo = object$fit$convInfo, control = object$fit$control, 
        na.action = object$fit$na.action, coefficients = param, parameters = param)
    if (correlation && rdf > 0) {
        ans$correlation <- cov2cor(vcov(object))
    }
    class(ans) <- c("summary.svynls","summary.nls")
    ans
}
