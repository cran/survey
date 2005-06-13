
calibrate<-function(design, ...) UseMethod("calibrate")

is.calibrated<-function(design){ !is.null(design$postStrata)}


calibrate.survey.design2<-function(design, formula, population,
                                   stage=NULL,  lambda=NULL,...){
  
  if (is.null(stage))
    stage<-if (is.list(population)) 1 else 0
  
  if(stage==0){
    ## calibration to population totals
    mm<-model.matrix(formula, model.frame(formula, design$variables))
    whalf<-sqrt(weights(design))
    sample.total<-colSums(mm*whalf*whalf)

    if (is.null(lambda))
      sigma2<-rep(1,nrow(mm))
    else
      sigma2<-drop(mm%*%lambda)
    
    if (length(sample.total)!=length(population))
      stop("Population and sample totals are not the same length.")
    
    if (any(names(sample.total)!=names(population)))
      warning("Sample and population totals have different names.")

    tqr<-qr(mm*whalf/sqrt(sigma2))
    if (is.null(lambda) && !all(abs(qr.resid(tqr,sigma2)) <1e-7))
      stop("Calibration models with constant variance must have an intercept")
    
    Tmat<-crossprod(mm*whalf/sqrt(sigma2))
    
    tT<-solve(Tmat,population-sample.total)
    
    g<-drop(1+mm%*%tT/sigma2)
    design$prob<-design$prob/g
    
    caldata<- list(qr=tqr, w=g*whalf*sqrt(sigma2), stage=0, index=NULL)
   
  } else {
    ## Calibration within clusters (Sarndal's Case C)
    if (stage>NCOL(design$cluster))
      stop("This design does not have stage",stage)

    if (!all(length(population[[1]])==sapply(population,length)))
      stop("Population totals are not all the same length")
    
    clusters<-unique(design$cluster[,stage])
    nc<-length(clusters)
    
    caldata<-list(qr=vector("list",nc), w=vector("list",nc),
                  stage=stage,index=as.character(clusters))

    mm<-model.matrix(formula, model.frame(formula, design$variables))

    if (is.null(lambda))
      sigma2<-rep(1,nrow(mm))
    else
      sigma2<-drop(mm%*%lambda)
    
    if(NCOL(mm)!=length(population[[1]]))
        stop("Population and sample totals are not the same length.")
      
    if (any(colnames(mm)!=names(population[[1]])))
      warning("Sample and population totals have different names.")
  
    stageweights<-1/apply(design$allprob[,1:stage,drop=FALSE],1,prod)
    if (any(duplicated(design$cluster[!duplicated(stageweights),stage])))
      stop("Weights at stage", stage, "vary within sampling units")

    cwhalf<-sqrt(weights(design)/stageweights)
    dwhalf<-sqrt(weights(design))
    tqr<-qr(mm)
    if (is.null(lambda) && !all(abs(qr.resid(tqr,sigma2)) <1e-7))
      stop("Calibration models with constant variance must have an intercept")
 
    for (i in 1:length(clusters)){ 
      cluster<-clusters[[i]]
      these<-which(cluster ==  as.character(design$cluster[,stage]))
      sample.total<-colSums(mm[these,,drop=FALSE]*cwhalf[these]*cwhalf[these])
      tqr<-qr(mm[these,,drop=FALSE]*cwhalf[these]/sqrt(sigma2[these]))
      Tmat<-crossprod(mm[these,,drop=FALSE]*cwhalf[these]/sqrt(sigma2[these]))
      tT<-solve(Tmat,population[[i]]-sample.total)
      g<-drop(1+mm[these,,drop=FALSE]%*%tT/sigma2[these])
      design$prob[these]<-design$prob[these]/g
      caldata$qr[[i]]<-tqr
      caldata$w[[i]]<-g*stageweights[these]*sqrt(sigma2[these])*cwhalf[these]^2
    }
  }  
  class(caldata)<-"greg_calibration"
  
  design$postStrata<-c(design$postStrata, list(caldata))
  design$call<-sys.call(-1)
  
  design
}


calibrate.svyrep.design<-function(design, formula, population,compress=NA,lambda=NULL,...){
  mf<-model.frame(formula, design$variables)
  mm<-model.matrix(formula, mf)
  whalf<-sqrt(design$pweights)

  if (is.null(lambda))
    sigma2<-rep(1,nrow(mm))
  else
    sigma2<-drop(mm%*%lambda)
  
  repwt<-as.matrix(design$repweights)
  if (!design$combined.weights)
    repwt<-repwt*design$pweights
  
  sample.total<-colSums(mm*whalf*whalf)

  if (length(sample.total)!=length(population))
    stop("Population and sample totals are not the same length.")
  if (any(names(sample.total)!=names(population)))
    warning("Sample and population totals have different names.")
  
  Tmat<-crossprod(mm*whalf)
  
  tT<-solve(Tmat,population-sample.total)
  
  g<-drop(1+mm%*%tT/sigma2)
  design$pweights<-design$pweights*g
  
  for(i in 1:NCOL(repwt)){
    whalf<-sqrt(repwt[,i])
    Tmat<-crossprod(mm*whalf/sqrt(sigma2))
    sample.total<-colSums(mm*whalf*whalf)
    g<-drop(1+mm%*%solve(Tmat,population-sample.total)/sigma2)
    repwt[,i]<-repwt[,i]*g
  }

  if (!design$combined.weights)
    repwt<-repwt/design$pweights

  if (compress ||
      (is.na(compress && inherits(design$repweights,"repweights_compressed")))){
    repwt<-compressWeights(repwt)
  }
    
  design$repweights<-repwt
  design$call<-sys.call(-1)

  design
}
