
calibrate<-function(design, ...) UseMethod("calibrate")

calibrate.survey.design2<-function(design, formula, population,
                                    aggregate.stage=NULL, stage=0, variance=NULL,
                                    bounds=c(-Inf,Inf), calfun=c("linear","raking","logit"),
                                    maxit=30, epsilon=1e-7, verbose=FALSE,
                                    ...){

  if (calfun=="linear" && (bounds==c(-Inf,Inf))){
    ## old code is better for ill-conditioned linear calibration
    rval<-regcalibrate(design,formula,population,
                       aggregate.stage=aggregate.stage, stage=stage,
                       lambda=variance,...)
    rval$call<-sys.call(-1)
    return(rval)
  }

  if (!is.null(aggregate.stage)){
    aggindex<-design$cluster[[aggregate.stage]]
  }

  expit<-function(x) 1-1/(1+exp(x))
  
  ## calibration to population totals
  mm<-model.matrix(formula, model.frame(formula, design$variables))
  ww<-weights(design)
  
  if (!is.null(aggregate.stage)){
    mm<-apply(mm,2,function(mx) ave(mx,aggindex))
    ww<-ave(ww,aggindex)
  }
  whalf<-sqrt(ww)
  sample.total<-colSums(mm*ww)
  
  if (length(sample.total)!=length(population))
    stop("Population and sample totals are not the same length.")
  
  if (any(names(sample.total)!=names(population)))
    warning("Sample and population totals have different names.")
  
  tqr<-qr(mm*whalf)
  if (is.null(!all(abs(qr.resid(tqr,whalf)))))
    warning("G-calibration models must have an intercept")

  calfun<-match.arg(calfun)
  g<-grake(mm,ww,calfun, bounds=bounds,population=population,
           verbose=verbose,epsilon=epsilon,maxit=maxit)
  
  design$prob<-design$prob/g
  
  caldata <- list(qr=tqr, w=g*whalf, stage=0, index=NULL)
  
  class(caldata) <- c("greg_calibration","gen_raking")
  
  design$postStrata <- c(design$postStrata, list(caldata))
  design$call <- sys.call(-1)
  
  design
}


calibrate.svyrep.design<-function(design, formula, population,compress=NA,
                                   aggregate.index=NULL, variance=NULL,
                                   bounds=c(-Inf,Inf), calfun=c("linear","raking","logit"),
                                   maxit=30, epsilon=1e-7, verbose=FALSE,
                                   ...){
  if (calfun=="linear" && (bounds==c(-Inf,Inf))){
    ## old code is better for ill-conditioned linear calibration
    rval<-regcalibrate(design,formula,population, compress=compress,
                       aggregate.index=aggregate.index,
                       lambda=variance,...)
    rval$call<-sys.call(-1)
    return(rval)
  }
  
  mf<-model.frame(formula, design$variables)
  mm<-model.matrix(formula, mf)
  ww<-design$pweights
  
  repwt<-as.matrix(design$repweights)
  if (!design$combined.weights)
    repwt<-repwt*design$pweights

  if (inherits(aggregate.index,"formula")){
    if (length(aggregate.index)!=2)
      stop("aggregate.index must be a one-sided formula")
    aggregate.index<-model.frame(aggregate.index, design$variables)
    if (NCOL(aggregate.index)>1)
      stop("aggregate.index must specify a single variable")
    aggregate.index<-aggregate.index[[1]]
  }
  
  if (!is.null(aggregate.index)){
    if (sqrt(max(ave(ww,aggregate.index,FUN=var),na.rm=TRUE))>1e-2*mean(ww))
      warning("Sampling weights are not constant within clusters defined by aggregate.index")
    mm<-apply(mm,2,function(mx) ave(mx,aggregate.index))
    ww<-ave(ww,aggregate.index)
    repwt<-apply(repwt,2,function(wx) ave(wx, aggregate.index))
  }
  whalf<-sqrt(ww)
  
  sample.total<-colSums(mm*ww)

  if (length(sample.total)!=length(population))
    stop("Population and sample totals are not the same length.")
  if (any(names(sample.total)!=names(population)))
    warning("Sample and population totals have different names.")
  

  calfun <- match.arg(calfun)
  gtotal <- grake(mm,ww,calfun,bounds=bounds,population=population,
                  verbose=verbose, epsilon=epsilon, maxit=maxit)

  design$pweights<-design$pweights*gtotal
  
  for(i in 1:NCOL(repwt)){
    wwi<-repwt[,i]
    if(verbose) cat("replicate = ",i,"\n")
    g<-grake(mm, wwi, calfun, eta=rep(0,NCOL(mm)), bounds=bounds, population=population,
             epsilon=epsilon, verbose=verbose, maxit=maxit)
    repwt[,i]<-as.vector(design$repweights[,i])*g
  }

  if (!design$combined.weights)
      repwt<-repwt/gtotal

  if (compress ||
      (is.na(compress && inherits(design$repweights,"repweights_compressed")))){
    repwt<-compressWeights(repwt)
  }
    
  design$repweights<-repwt
  design$call<-sys.call(-1)

  design
}



grake<-function(mm,ww,calfun,eta=rep(0,NCOL(mm)),bounds,population,epsilon, verbose,maxit){

  sample.total<-colSums(mm*ww)
  require(MASS) ##ginv
  
  Fm1<-switch(calfun,
              linear=function(u) pmin(pmax(u+1,bounds[1]),bounds[2])-1,
              raking=function(u) pmin(pmax(exp(u),bounds[1]),bounds[2])-1,
              logit=function(u) {
                L <- bounds[1]
                U <- bounds[2]
                A <- (U-L)/((U-1)*(1-L))
                eAu <- exp(A*u)
                ( L*(U-1) + U*(1-L)*eAu)/(U-1+(1-L)*eAu)-1
              } 
              )
  dF<-switch(calfun,
             linear=function(u) as.numeric(u<bounds[2]-1 & u>bounds[1]-1),
             raking=function(u) ifelse(u<bounds[2]-1 & u>bounds[1]-1,exp(u),0),
             logit=function(u) {
               L <- bounds[1]
               U <- bounds[2]
               A <- (U-L)/((U-1)*(1-L))
               eAu <- exp(A*u)
               U*(1-L)*eAu*A/(U-1+(1-L)*eAu)-(L*(U-1)+U*(1-L)*eAu)/(U-1+(1-L)*eAu)^2
             }
             )

  xeta<-drop(mm%*%eta)
  g<-1+Fm1(xeta)

  iter<-1

  repeat({
    Tmat<-crossprod(mm*ww*dF(xeta), mm)

    misfit<-(population-sample.total-colSums(mm*ww*Fm1(xeta)))
    deta<-ginv(Tmat)%*%misfit
    eta<-eta+deta

    xeta<- drop(mm%*%eta)
    g<-1+Fm1(xeta)
    misfit<-(population-sample.total-colSums(mm*ww*Fm1(xeta)))
    
    if (verbose)
      print(misfit)

    if (all(abs(misfit)/(1+abs(population))<epsilon)) break

    iter <- iter+1
    if (iter>maxit) stop("Failed to converge in ",iter," iterations")
  })

  attr(g,"eta")<-eta
  g
}
