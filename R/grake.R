
make.calfun<-function(Fm1,dF, name){
  if (!identical(names(formals(Fm1)), c("u","bounds")))
    stop("wrong argument names for Fm1")
  if(!identical(names(formals(dF)), c("u","bounds")))
    stop("wrong argument names for dF")
  rval<-list(Fm1=Fm1, dF=dF, name=name)
  class(rval)<-"calfun"
  rval
}


cal_names <- function(formula, design, ...) UseMethod("cal_names",design)

cal_names.DBIsvydesign<-function(formula, design,...){
    design$variables <- getvars(formula, design$db$connection, 
        design$db$tablename, updates = design$updates, subset = design$subset)
    colnames(model.matrix(formula, model.frame(formula,design$variables[0,])))
}

cal_names.ODBCsvydesign<-function(formula, design,...){
    design$variables <- getvars(formula, design$db$connection, 
                                design$db$tablename, updates = design$updates)
    colnames(model.matrix(formula, model.frame(formula,design$variables[0,])))
}

cal_names.survey.design<-function(formula,design,...) {
    colnames(model.matrix(formula, model.frame(formula,model.frame(design)[0,])))
}
    
print.calfun<-function(x,...) cat("calibration metric: ",x$name,"\n")

calibrate<-function(design, ...) UseMethod("calibrate")

calibrate.survey.design2<-function(design, formula, population,
                                    aggregate.stage=NULL, stage=0, variance=NULL,
                                    bounds=c(-Inf,Inf), calfun=c("linear","raking","logit"),
                                    maxit=50, epsilon=1e-7, verbose=FALSE, force=FALSE, trim=NULL,
                                    bounds.const = FALSE, sparse=FALSE,
                                    ...){

  if(is.list(formula) && is.list(population)){
    ## inputs as marginal totals, as in rake()
    population<-margins2totals(formula,population)
    formula<-as.formula(paste("~",paste(sapply(formula,function(f) paste(all.vars(f),collapse="*")),collapse="+")))
    if (verbose){
      print(formula)
      print(population)
    }
  }

  
  if (!is.list(bounds)) bounds <- list(lower = bounds[1], upper = bounds[2]) ## converting to list if not already
  if (is.character(calfun)) calfun<-match.arg(calfun)
  if (is.character(calfun) && calfun=="linear" && all(unlist(bounds) == c(-Inf, Inf))){
    ## old code is better for ill-conditioned linear calibration
    rval<-regcalibrate(design,formula,population,
                       aggregate.stage=aggregate.stage, stage=stage,
                       lambda=variance,sparse=sparse,...)
    rval$call<-sys.call(-1)
    return(rval)
  }

  if(is.character(calfun))
    calfun<-switch(calfun,linear=cal.linear, raking=cal.raking, logit=cal.logit)
  else
    if(!inherits(calfun,"calfun"))
      stop("'calfun' must be a string or of class 'calfun'.")

  if (length(epsilon)!=1 && length(epsilon)!=length(population))
    stop("'epsilon' must be a scalar or of the same length as 'population'")
  
  if (!is.null(aggregate.stage)){
    aggindex<-design$cluster[[aggregate.stage]]
  }

  expit<-function(x) 1-1/(1+exp(x))
  
  ## calibration to population totals
  if(sparse){
    mm<-sparse.model.matrix(formula, model.frame(formula, model.frame(design)))
  }else{
    mm<-model.matrix(formula, model.frame(formula, model.frame(design)))
  }
  
  ww<-weights(design)
  
  if (bounds.const) bounds<-lapply(bounds,
        function(x) ifelse(ww>0,x/ww,x)) # if bounds are set to a constant convert to multiplicative value and run as normal
  
  if (!is.null(aggregate.stage)){
    mm<-apply(mm,2,function(mx) ave(mx,aggindex))
    ww<-ave(ww,aggindex)
  }
  whalf<-sqrt(ww)
  sample.total<-colSums(mm*ww)
  
  if(any(sample.total==0)){
    ## drop columsn where all sample and population are zero
    zz<-(population==0) & (apply(mm,2,function(x) all(x==0)))
    mm<-mm[,!zz]
    population<-population[!zz]
    sample.total<-sample.total[!zz]
    if (length(epsilon)>1) epsilon <- epsilon[!zz]
  }

    
  if (length(sample.total)!=length(population)){
    print(sample.total)
    print(population)
    stop("Population and sample totals are not the same length.")
  }
  if(!is.null(names(population))){
    if (!all(names(sample.total) %in% names(population))){
        warning("Sampling and population totals have different names.")
        cat("Sample: "); print(names(sample.total))
        cat("Popltn: "); print(names(population))
        }
    else if (!all(names(sample.total) == names(population))){
      warning("Sample and population totals reordered to make names agree: check results.")
      population <- population[match(names(sample.total), names(population))]
    }
  }
  
  tqr<-qr(mm*whalf)
  #if (!all(abs(qr.resid(tqr,whalf))<1e-10))
  #  warning("G-calibration models must have an intercept")

  g<-grake(mm,ww,calfun, bounds=bounds,population=population,
           verbose=verbose,epsilon=epsilon,maxit=maxit, variance=variance)

  if(!is.null(trim)) {
      gnew<-pmax(trim[1], pmin(g, trim[2]))
      outside<-g<trim[1] | g>trim[2]
      if (any(outside)){
        trimmings<-(g-gnew)*ww
        gnew[!outside]<-gnew[!outside]+sum(trimmings)/sum(ww[!outside])
        g<-gnew
        attr(g,"failed")<-NULL
        message(paste(sum(outside),"weights were trimmed"))
      }
    }
  if (!is.null(attr(g,"failed"))){
    if (!force) stop("Calibration failed")
  }
  
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
                                   maxit=50, epsilon=1e-7, verbose=FALSE,force=FALSE, trim=NULL,
                                   bounds.const=FALSE, sparse=FALSE,
                                   ...){

  if(is.list(formula) && is.list(population)){
    ## inputs as marginal totals, as in rake()
    population<-margins2totals(formula,population)
    formula<-as.formula(paste("~",paste(sapply(formula,function(f) paste(all.vars(f),collapse="*")),collapse="+")))
    if (verbose){
      print(formula)
      print(population)
    }
  }
   
  if (is.character(calfun)) calfun<-match.arg(calfun)
  if (length(epsilon)!=1 && length(epsilon)!=length(population))
    stop("'epsilon' must be a scalar or of the same length as 'population'")
  
  if (!is.list(bounds)) bounds <- list(lower = bounds[1], upper = bounds[2]) ## converting to list if not already
  if (is.character(calfun) && calfun=="linear" && all(unlist(bounds)==c(-Inf,Inf))){
    ## old code is better for ill-conditioned linear calibration
    rval<-regcalibrate(design,formula,population, compress=compress,
                       aggregate.index=aggregate.index,
                       lambda=variance,...)
    rval$call<-sys.call(-1)
    return(rval)
  }
  
  mf<-model.frame(formula, design$variables)
  if(sparse){
    mm<-sparse.model.matrix(formula, model.frame(formula, model.frame(design)))
  }else{
    mm<-model.matrix(formula, model.frame(formula, model.frame(design)))
  }
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
  
  if(any(sample.total==0)){
    ## drop columsn where all sample and population are zero
    zz<-(population==0) & (apply(mm,2,function(x) all(x==0)))
    mm<-mm[,!zz]
    population<-population[!zz]
    sample.total<-sample.total[!zz]
    if (length(epsilon)>1) epsilon <- epsilon[!zz]
  }
  
    if (length(sample.total)!=length(population)){
        print(sample.total)
        print(population)
        stop("Population and sample totals are not the same length.")
    }
    if (!is.null(names(population))){
        if (!all(names(sample.total) %in% names(population))){
            warning("Sampling and population totals have different names.")
            cat("Sample: "); print(names(sample.total))
            cat("Popltn: "); print(names(population))
        }
        else if (!all(names(sample.total) == names(population))){
            warning("Sample and population totals reordered to make names agree: check results.")
            population <- population[match(names(sample.total), names(population))]
        }
    }
    
  if(is.character(calfun))
    calfun<-switch(calfun, linear=cal.linear, raking=cal.raking, logit=cal.logit)
  else if (!inherits(calfun,"calfun"))
    stop("'calfun' must be a string or a 'calfun' object")
  gtotal <- grake(mm,ww,calfun,bounds=bounds,population=population,
                  verbose=verbose, epsilon=epsilon, maxit=maxit, 
                  variance=variance)

  if(!is.null(trim)) {
      gnew<-pmax(trim[1], pmin(gtotal, trim[2]))
      outside<-gtotal<trim[1] | gtotal>trim[2]
      if (any(outside)){
        trimmings<-(gtotal-gnew)*ww
        gnew[!outside]<-gnew[!outside]+sum(trimmings)/sum(ww[!outside])
        gtotal<-gnew
        attr(gtotal,"failed")<-NULL
        message(paste(sum(outside),"weights were trimmed"))
      }
    }
  if (!force && !is.null(attr(gtotal,"failed"))) stop("Calibration failed")

  design$pweights<-design$pweights*gtotal
  
  for(i in 1:NCOL(repwt)){
    wwi<-repwt[,i]
    if(verbose) cat("replicate = ",i,"\n")
    g<-grake(mm, wwi, calfun, eta=rep(0,NCOL(mm)), bounds=bounds, population=population,
             epsilon=epsilon, verbose=verbose, maxit=maxit, 
             variance=variance)
    
    if(length(trim)==2){
      outside<-(g<trim[1]) | (g>trim[2])
      if (any(outside)) {
        gnew<-pmax(trim[1],pmin(g,trim[2]))
        trimmings<-(g-gnew)*wwi
        gnew[!outside]<-gnew[!outside]+sum(trimmings)/sum(wwi[!outside])
        g<-gnew
      }}

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

cal.linear<-make.calfun(function(u,bounds) pmin(pmax(u+1,bounds$lower),bounds$upper)-1,
                        function(u, bounds) as.numeric(u < bounds$upper-1 & u > bounds$lower-1),
                        "linear calibration")
cal.raking<-make.calfun(function(u,bounds) pmin(pmax(exp(u),bounds$lower),bounds$upper)-1,
                        function(u, bounds) ifelse(u<bounds$upper-1 & u>bounds$lower-1,exp(u),0),
                        "raking")
cal.logit<-make.calfun(
  function(u,bounds) {
    if (any(!is.finite(unlist(bounds)))) stop("Logit calibration requires finite bounds")
    L <- bounds$lower
    U <- bounds$upper
    A <- (U-L)/((U-1)*(1-L))
    eAu <- exp(A*u)
    ( L*(U-1) + U*(1-L)*eAu)/(U-1+(1-L)*eAu)-1
  },
  function(u,bounds) {
    L <- bounds$lower
    U <- bounds$upper
    A <- (U-L)/((U-1)*(1-L))
    eAu <- exp(A*u)
    U*(1-L)*eAu*A/(U-1+(1-L)*eAu)-( (L*(U-1)+U*(1-L)*eAu)*( (1-L)*eAu*A ) )/(U-1+(1-L)*eAu)^2
  },
  "logit calibration"
)

cal.sinh <- make.calfun(
  Fm1 = function(u, bounds) 
  {
    if (any(!is.finite(unlist(bounds)))) 
      stop("Sinh calibration requires finite bounds")

    
    dasinh <- function(u, alpha = 1) {
      p <- asinh(2 * alpha * u)
      m <- (p + sqrt((p^2) + 4*(alpha^2)))/(2*alpha)
      m <- as.vector(m)
      m
    }
    # linear truncation
    pmin(pmax(dasinh(u), bounds$lower), bounds$upper) - 1
    
  },
  dF = function(u, bounds)
  {
    ddasinh <- function(u, alpha = 1) {
      p <- asinh(2 * alpha * u)
      m <- 1/sqrt(1 + (2*alpha*u)^2)*(1 + p/sqrt(p^2 + 4*alpha^2))
      m <- as.vector(m)
      m
    }
    
    ifelse(u < bounds$upper - 1 & u > bounds$lower - 1, ddasinh(u), 0)
    
  },
  name = 'sinh calibration'
  
)
                      
grake<-function(mm,ww,calfun,eta=rep(0,NCOL(mm)),bounds,population,epsilon, 
                verbose, maxit, variance=NULL){
  sample.total<-colSums(mm*ww)
  if(!inherits(calfun,"calfun")) stop("'calfun' must be of class 'calfun'")
  
  Fm1<-calfun$Fm1
  dF<-calfun$dF
  
  if (is.null(variance)){
    sigma2<-rep(1,nrow(mm))
  }
  else if(length(variance) == nrow(mm)){
    sigma2<-drop(variance)
  }else{
    sigma2<-drop(mm%*%variance)
  }
  
  xeta<-drop(mm%*%eta/sigma2)
  g<-1+Fm1(xeta, bounds)
  deriv <- dF(xeta, bounds)
  iter<-1

  ## pre-scaling for people starting with no weights
  SOMETHRESHOLD<-20
  scales<-population/sample.total
  if (min(scales)> SOMETHRESHOLD){
    scale<-mean(scales)
    ww<-ww*scale
    sample.total<-sample.total*scale
    if(verbose) message(paste("Sampling weights rescaled by",signif(scale,3)))
    if (any(is.finite(unlist(bounds)))) warning(paste("Bounds were set but will be interpreted after rescaling by",signif(scale,3)))
  } else scale<-NULL
  
  repeat({
    Tmat<-crossprod(mm*ww/sqrt(sigma2)*deriv, mm/sqrt(sigma2))

    misfit<-(population-sample.total-colSums(mm*ww*Fm1(xeta, bounds)))
    deta<-MASS::ginv(as(Tmat, "matrix"), tol=256*.Machine$double.eps)%*%misfit
    eta<-eta+deta

    xeta<- drop(mm%*%eta/sigma2)
    g<-1+Fm1(xeta, bounds)
    deriv <- dF(xeta, bounds)
    while(iter<maxit && any(!is.finite(g),!is.finite(deriv))){
      iter<-iter+1
      deta<-deta/2
      eta<-eta-deta
      xeta<- drop(mm%*%eta/sigma2)
      g<-1+Fm1(xeta, bounds)
      deriv <- dF(xeta, bounds)
      if(verbose) print("Step halving")
    }
    misfit<-(population-sample.total-colSums(mm*ww*Fm1(xeta, bounds)))
    
    if (verbose)
      print(misfit)

    if (all(abs(misfit)/(1+abs(population))<epsilon)) break

    iter <- iter+1
    if (iter>maxit) {
       achieved<-max((abs(misfit)/(1+abs(population))))
       warning("Failed to converge: eps=",achieved," in ",iter," iterations")
       attr(g,"failed")<-achieved
       break;
     }
  })

  if (!is.null(scale)) g<-g*scale
  attr(g,"eta")<-eta
  g
}


trimWeights<-function(design, upper=Inf,lower=-Inf, ...){
  UseMethod("trimWeights")
}


trimWeights.survey.design2<-function(design, upper=Inf, lower= -Inf, strict=FALSE,...){
  pw<-weights(design,"sampling")
  outside<-pw<lower | pw>upper
  if (!any(outside)) return(design)
  pwnew<-pmax(lower,pmin(pw, upper))
  trimmings<-pw-pwnew
  pwnew[!outside]<-pwnew[!outside]+sum(trimmings)/sum(!outside)
  design$prob<-1/pwnew
  design$call<-sys.call()
  design$call[[1]]<-as.name(.Generic)
  if (strict) ## ensure that the trimmings don't push anything outside the limits
    trimWeights(design, upper,lower, strict=TRUE)
  else 
    design
}

trimWeights.svyrep.design<-function(design, upper=Inf, lower= -Inf, compress=FALSE,...){
  pw<-weights(design,"sampling")
  outside<-pw<lower | pw>upper
  if (any(outside)) {
    pwnew<-pmax(lower,pmin(pw, upper))
    trimmings<-pw-pwnew
    pwnew[!outside]<-pwnew[!outside]+sum(trimmings)/sum(!outside)
    design$prob<-1/pw
  }
  rw<-weights(design, "analysis")
  outside<-rw<lower | rw>upper
  if (any(outside)) {
    rwnew<-pmax(lower,pmin(rw, upper))
    trimmings<-rw-rwnew
    rwnew<-rwnew[!outside]+t(t(!outside)+colSums(trimmings)/colSums(!outside))
    if (compress)
      design$repweights<-compressWeights(rwnew)
    else
      design$repweights<-rwnew
    design$combined.weights<-TRUE
  }
  
  design
}


margins2totals<-function(formulas, totals){
	totals<-mapply(onemargin2totals,formulas,totals,SIMPLIFY=FALSE)
	totaln<-do.call(c,totals)
	totalorder<-do.call(c,lapply(totals,function(x) attr(x,"order")))
	totaln<-totaln[order(totalorder)]

	totaln[!duplicated(names(totaln))]	
	}
	
	
onemargin2totals<-function(formula,total){
    if (is.table(total)) total<-as.data.frame(total)
    if (!is.data.frame(total) && is.vector(total) && (length(formula[[2]])==1)){
        ## just a vector
        total<-as.table(total)
        d<-dimnames(total)
        names(d)<-paste(deparse(formula[[2]]),collapse="")
        total<-as.data.frame(total)
    }
    if (!is.data.frame(total)) stop("incorrect format for population totals")
    
    newformula<-as.formula(paste("Freq",paste(all.vars(formula),collapse="*"),sep="~"))
    mf<-model.frame(newformula,as.data.frame(total))
    mm<-model.matrix(newformula,mf)
    intorder<-c(1,attr(terms(newformula),"order")[attr(mm,"assign")])
    rval<-colSums(mf$Freq*mm)
    attr(rval,"order")<-intorder
    rval
}	
