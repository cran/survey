
reweight<-function(design,...) UseMethod("reweight")

reweight.dualframe<-function(design, targets=NULL, totals=NULL,
                             estimator=c("constant","expected"), theta=NULL,
                             theta_grid=seq(0,1,by=0.05),...) {
    
    estimator<-match.arg(estimator)
    nframes<-length(design$designs)
    overlaps<-design$overlaps
    frame_weights<-vector("list",nframes)
    
    if (estimator=="expected"){
        maxes<-sapply(overlaps, max, na.rm=TRUE)
        mins<-sapply(overlaps, function(x) min(x[x>0], na.rm=TRUE))
        if (max(mins)>=1) 
            overlaps_are_weights<-TRUE
        else if(min(maxes)>1)
            stop("overlaps must be all >=1 or all <=1")
        else
            overlaps_are_weights<-FALSE
        
        for(f in 1:nframes){
            overlaps[[f]][overlaps[[f]]==0]<-NA
            if (overlaps_are_weights) overlaps[[f]]<-1/overlaps[[f]]
            frame_weights[[f]]<-(1/rowSums(1/overlaps[[f]],na.rm=TRUE))/weights(design$designs[[f]],"sampling")
        }
        frame_weights<-do.call(c, frame_weights)
        
        design_weights<-do.call(c, lapply(design$designs, weights, type="sampling"))
        
        rval<-list(designs=design$designs,overlaps=design$overlaps, frame_scale=design$frame_scale,
                   frame_weights=frame_weights, design_weights=design_weights,
                   call=sys.call(), dchecks=design$dchecks, estimator=estimator)
        class(rval)<-class(design)
        return(rval)  ##done
    }
    if (estimator=="constant" && !is.null(theta)){
        frame_scale<-c(theta, 1-theta)
        for(f in 1:nframes){
            frame_weights[[f]]<-ifelse(overlaps[[f]][,3-f]>0, frame_scale[f], 1)
        }
        frame_weights<-do.call(c, frame_weights)
        design_weights<-do.call(c, lapply(design$designs, weights, type="sampling"))
        
        rval<-list(designs=design$designs,overlaps=design$overlaps, frame_scale=design$frame_scale,
                   frame_weights=frame_weights,  design_weights=design_weights,
                   call=sys.call(), dchecks=design$dchecks, estimator=estimator)
        class(rval)<-class(design)
        return(rval) ##done
    }
    
    ## if we get here we are optimising something
    if (!xor(is.null(targets), is.null(totals)))
        stop("for estimator='constant', must provide exactly one of targets, totals, and theta")

    if (!is.null(totals)){
        targets<-lapply(totals,
                        function(formula) bquote(vcov(svytotal(.(formula), design=.DESIGN,na.rm=TRUE))))
        }
    if(!is.null(design$theta))
        theta_old<-design$theta
    else
        theta_old<-1/length(design$designs)
    
    ntargets<-length(targets)
    nthetas<-length(theta_grid)
    variances<-lapply(1:ntargets, function(x) numeric(nthetas))   
    for (j in 1:nthetas){
        frame_scale<-c(theta_grid[j], 1-theta_grid[j])
        for(i in 1:ntargets){
            frame_weights<-vector("list",nframes)
            for(f in 1:nframes){
                frame_weights[[f]]<-ifelse(overlaps[[f]][,3-f]>0, frame_scale[f], 1)
            }
            frame_weights<-do.call(c, frame_weights)
            design_weights<-do.call(c, lapply(design$designs, weights, type="sampling"))
            
            tempval<-list(designs=design$designs,overlaps=overlaps, frame_scale=design$frame_scale,
                          frame_weights=frame_weights,  design_weights=design_weights,
                          call=sys.call(), dchecks=design$dchecks, estimator=estimator)
            class(tempval)<-class(design)
            target<-do.call(substitute, list(targets[[i]], list(.DESIGN=tempval)))
            estimator<-eval(target)
            if (length(estimator)>1) {
                warning(paste("multiple variances reported, only the first used for",deparse(targets[[i]])))
                estimator<-estimator[1]
                }
            variances[[i]][j]<-estimator
        }
    }
    reweight_info<-list(theta=theta_grid, targets=targets, variances=variances,
                        theta_old=theta_old, opt_thetas=theta_grid[sapply(variances, which.min)])
    class(reweight_info)<-"reweight_info"
    rval<-design
    rval$rewt<-reweight_info
    class(rval)<-c("dualframe_with_rewt", class(design))
    rval
}

coef.dualframe_with_rewt<-function(object, ...) object$rewt$opt_thetas

plot.dualframe_with_rewt<-function(x,y,type="b",...){
    ntargets<-length(x$rewt$targets)
    scaled_vars<-vector("list",ntargets)
    
    if((x$rewt$theta_old < min(x$rewt$theta)) || (x$rewt$theta_old > max(x$rewt$theta)))
        theta_old<-median(x$rewt$theta)
    else
        theta_old<-x$rewt$theta_old
    
    for(i in 1:ntargets){
        fn<-approxfun(x$rewt$theta, x$rewt$variances[[i]])
        old_var<-fn(theta_old)
        scaled_vars[[i]]<-x$rewt$variances[[i]]/old_var
    }
    matplot(x$rewt$theta, do.call(cbind,scaled_vars),type=type,xlab="theta",ylab="Scaled variance",...)
    
    invisible(list(x$rewt$theta,scaled_vars))
}

reweight.multiframe<-function(design, ...) stop("multi-frame reweighting under construction")

