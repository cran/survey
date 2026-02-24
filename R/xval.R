

withCrossval <- function(design, formula, trainfun, testfun,
                         loss=c("MSE","entropy","AbsError"),
                         intercept,tuning=NULL, nearly_zero=1e-4,...){
    UseMethod("withCrossval")
}

## we may want to cluster the repweights to keep the number down. 
##
withCrossval.svyrep.design<-function(design, formula, trainfun, testfun,
                                     loss=c("MSE","entropy","AbsError"),
                                     intercept, tuning, nearly_zero=1e-4,...){

    if (is.character(loss)) {
        loss<-match.arg(loss)
        lossfun<-switch(loss, MSE=function(y,hat,w) sum(w*(y-hat)^2),
                        AbsError=function(y,hat,w) sum(w*abs(y-hat)),
                        entropy=function(y, hat, w) sum(w*(y*log(hat)+(1-y)*log(1-hat))))
    } else if (!is.function(loss)){
        lossfun<-loss
    }
       
    repweights<-weights(design, "analysis")
    pweights<-weights(design,"sampling")

    testset<- (repweights/pweights)<= nearly_zero
    if (any(colSums(testset)==0))
        stop("some replicates have no test-set observations. Check nearly_zero.")
    
    mf<-model.frame(formula, model.frame(design))
    y<-model.response(mf)
    X<-model.matrix(formula, mf)
        
    hat<-matrix(NA, ncol=ncol(repweights),nrow=nrow(repweights))
    losses<-numeric(length(tuning))

    for(i in 1:length(tuning)){
        for(fold in 1:ncol(repweights)){
            is_test<-testset[,fold]
            is_train<-!testset[,fold]
            w<-repweights[,fold]
            fit<-trainfun(X[is_train,,drop=FALSE], y[is_train], w[is_train], tuning=tuning[i])
            hat[is_test,fold]<-testfun(X[is_test,], trainfit=fit,tuning=tuning[i])
        }

        loss<-0
        for(fold in 1:ncol(hat)){
            noNA<-!is.na(hat[,fold])
            loss<-loss+lossfun(y[noNA],hat[noNA,fold],pweights[noNA])*design$rscales[fold]
        }
        losses[i]<-loss*design$scale
    }

    losses/sum(pweights)
}

