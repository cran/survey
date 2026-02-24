### multiphase

na_failsafe<-function(message="missing values in object"){
    function(object,...){
        if (NCOL(object)==0)
            object
        else {
            ok <- complete.cases(object)
            if (all(ok)) 
                object
            else stop(message)
        }
    }
}

## less memory-hungry version for sparse tables
svy_interaction<-function (..., drop = TRUE) {
    args <- list(...)
    narg <- length(args)
    if (narg == 1 && is.list(args[[1]])) {
        args <- args[[1]]
        narg <- length(args)
    }
    ls<-sapply(args,function(a) length(levels(a)))
    ans<-do.call("paste",c(lapply(args,as.character),sep="."))
    ans<-factor(ans)
    return(ans)
}


lookup_ids<-function(ids, subsets, data){
    na_id<-na_failsafe("missing values in `id'")
    if(is.null(ids))
        stop("must provide ids")

    m<-length(ids)
    if(m<=1 || !is.list(ids))
        stop("multiphase needs at least two phases")

    for(i in 1:m){
        istar<-max(i-1,1)
        if(inherits(ids[[i]],"formula")) {
            mf<-substitute(model.frame(ids[[i]], data=data[subsets$cumulative[[istar]],,drop=FALSE], na.action=na_id))
            ids[[i]]<-eval.parent(mf)
            if (ncol(ids[[i]])==0) ## formula was ~1
                ids[[i]]<-data.frame(id=1:nrow(ids[[i]])[subsets$cumulative[[istar]]])
        } else{
            if (is.null(ids[[i]]))
                stop("Must provide ids at each phase")
            else
                ids[[i]]<-na_id(data.frame(ids[[i]]))
        }
        
        ## make ids factor if they are character
        for(j in 1:ncol(ids[[i]])){
            if (is.character(ids[[i]][[j]]))
                ids[[i]][[j]]<-factor(ids[[i]][[j]])
        }
    }
    ids
    
}

na_subset<-function(message="missing values not allowed in subset", subset=NULL){
    function(object){
         if (NCOL(object)==0)
            object
         else {
             ok <- complete.cases(object)
             if (!is.null(subset))
                 ok<-ok[subset]
             if (all(ok)) 
                 object
             else stop(message)
         }

    }
}


lookup_probs<-function(probs,subsets, data){

    m<-length(probs)
    if (m<=1) stop("need probs (or NULL) at each phase")
    
    for (i in 1:m){
        this_subset<-subsets$cumulative[[i]]
        na.prob<-na_subset("missing values in `prob'", this_subset)
        if(is.null(probs[[i]])) next
        if(inherits(probs[[i]],"formula")){
            mf<-substitute(model.frame(probs[[i]],data=data,na.action=na.prob))
            probs[[i]]<-eval.parent(mf)
            if (ncol(probs[[i]])==0) probs[[i]]<-data.frame(prob=rep(1, sum(this_subset)))
        } else  if(!inherits(probs[[i]],"pps_spec")) stop("prob must be list of formulas or pps_spec objects")
    }
    probs   
}


compute_probs<-function(strata, cluster, subset){
        popsize <- mapply(function(s, i) ave(!duplicated(i), s, FUN = sum),
                          strata, cluster)
        
        fpc <- as.fpc(popsize[subset,,drop=FALSE], strata[subset,,drop=FALSE],cluster[subset,,drop=FALSE])
        allprob <- 1/weights(fpc, final = FALSE)
        prob <- apply(as.data.frame(allprob), 1, prod)
        list(allprob=allprob,prob=prob)     
    }



lookup_subset<-function(subset, data){
    ## Note: m is nphases-1 here: two phases, one subset
    m<-length(subset)
    cumulative<-vector("list",m+1)
    cumulative[[1]] <-rep(TRUE, NROW(data))
    
    for(i in 1:m){
        na.prob<-na_subset("missing values in `subset'", cumulative[[i]])

        if(inherits(subset[[i]],"formula")){
            mf<-substitute(model.frame(subset[[i]],data=data,na.action=na.prob)[[1]])
            subset[[i]]<-eval.parent(mf)
        } else  stop("subset must be a list of formulas ")
        added<-isTRUE(subset) & !cumulative[[i]]
        if(any(added)) stop(paste("subsets must be nested; failed for",subset))
        subset[[i]][is.na(subset[[i]])]<-FALSE
        cumulative[[i+1]]<-cumulative[[i]] & subset[[i]]
    }
    return(list(each=subset, cumulative=cumulative))
}

lookup_strata<-function(strata,ids, subsets, data){
    m<-length(strata)

    for(i in 1:m){
        istar<-max(i-1,1)
        na_strata<-na_subset("missing values in `strata'", subsets$cumulative[[istar]])

        if (!is.null(strata[[i]])){
            if(inherits(strata[[i]],"formula")){
                mf<-substitute(model.frame(strata[[i]], data=data[subsets$cumulative[[istar]],,drop=FALSE],
                                           na.action=na.fail))
                strata[[i]]<-eval.parent(mf)
            }
            if (!is.list(strata[[i]]))
                strata[[i]]<-data.frame(strata=strata[[i]])
            
            for(k in 1:NCOL(strata[[i]])){ ##drop empty strata
                if (is.factor(strata[[i]][[k]])) {
                    strata[[i]][[k]]<-as.factor(as.character(strata[[i]][[k]]))
                } else if (is.character(strata[[i]][[k]])){ ##coerce string to factor
                    strata[[i]][[k]]<-as.factor(strata[[i]][[k]])
                }
            }
        } else {
            strata[[i]]<-na_strata(as.data.frame(matrix(1, nrow=NROW(ids[[i]]), ncol=NCOL(ids[[i]]))))
        }
    }
    
    strata
}


check_multistage<-function(ids, strata){
    ## check multistage setup for one phase
    
    has.strata<-NROW(unique(strata))>1
    ## check for only one PSU: probably a typo
    if ((length(unique(ids[,1]))==1) && !has.strata){
        stop("Phase has only one primary sampling unit")
    }
    
      ## force subclusters nested in clusters
      if (NCOL(ids)>1){
        N<-ncol(ids)
        for(i in 2:N){
          ids[,i]<-do.call("svy_interaction", ids[,1:i,drop=FALSE])
        }
      }
      
    ## check if clusters nested in strata
    ## nest=TRUE not allowed here
    if(!is.null(strata) && NCOL(ids)){
       sc<-(rowSums(table(ids[,1],strata[,1])>0))
       if(any(sc>1)) stop("Clusters not nested in strata at top level")
    }

    ## force substrata nested in clusters
    N<-ncol(ids)
    NS<-ncol(strata)
    if (N>1){
        for(i in 2:N)
            strata[,i]<-interaction(strata[,min(i,NS)], ids[,i-1])
    }
    
    list(ids=ids, strata=strata, has.strata=has.strata)
}


phase_for_var<-function(data, subsets){
    ## yes, this could be optimised, but it runs once
    ## and nvars*nphases isn't going to be worse than thousands
    nvars<-ncol(data)
    variable_phases<-rep(1,nvars)
    names(variable_phases)<-names(data)
    
    nphases<-length(subsets)
    for(i in 2:nphases){
        this_phase<-subsets$cumulative[[i]]
        for(j in 1:nvars){
            if (all(is.na(data[!this_phase,j])))
                variable_phases[j]<-i
        }
    }
    variable_phases
}
 

weights.multiphase<-function(object, type=c("sampling","phase"),...){
    type<-match.arg(type)
    nphases<-object$nphases
    if(type=="sampling") {
        w<-object$finalweights*object$subpop[object$subsets$cumulative[[nphases]],nphases]
    } else if(type=="phase"){
        w<-vector("list",nphases)
        for(i in 1:nphases){
            w[[i]]<-object$phaseweights[[i]]*object$subpop[object$subsets$cumulative[[i]],i]
        }
      
    } else stop("can't happen")
    w
}


update.multiphase<-function(object,...){
    ## use all available data, so that rephase() can work
    dots <- substitute(list(...))[-1]
    newnames <- names(dots)
    for (j in seq(along = dots)) {
        object$data[, newnames[j]] <- eval(dots[[j]], object$variables, 
            parent.frame())
    }
    object$call <- sys.call(-1)
    object
}
    


multiphase_getdata<-function(x, design, formula_only=FALSE, phase=design$nphases){
    if (inherits(x, "formula")) {
        mf <- model.frame(x, design$data,  na.action = na.pass)
        xx <- lapply(attr(terms(x), "variables")[-1], function(tt) model.matrix(eval(bquote(~0 + 
            .(tt))), mf))
        cols <- sapply(xx, NCOL)
        x <- matrix(nrow = NROW(xx[[1]]), ncol = sum(cols))
        scols <- c(0, cumsum(cols))
        for (i in 1:length(xx)) {
            x[, scols[i] + 1:cols[i]] <- xx[[i]]
        }
        colnames(x) <- do.call("c", lapply(xx, colnames))
    }
    else {
        if (formula_only) stop("needs a model formula")
        if (typeof(x) %in% c("expression", "symbol")) 
            x <- eval(x, design$data)
        else {
            if (is.data.frame(x) && any(sapply(x, is.factor))) {
                xx <- lapply(x, function(xi) {
                  if (is.factor(xi)) 
                    0 + (outer(xi, levels(xi), "=="))
                  else xi
                })
                cols <- sapply(xx, NCOL)
                scols <- c(0, cumsum(cols))
                cn <- character(sum(cols))
                for (i in 1:length(xx)) cn[scols[i] + 1:cols[i]] <- paste(names(x)[i], 
                  levels(x[[i]]), sep = "")
                x <- matrix(nrow = NROW(xx[[1]]), ncol = sum(cols))
                for (i in 1:length(xx)) {
                  x[, scols[i] + 1:cols[i]] <- xx[[i]]
                }
                colnames(x) <- cn
            }
        }
        if(NROW(x)==sum(design$subsets$cumulative[[phase]]))
            return(as.matrix(x))
        if(NROW(x)==sum(design$subsets$cumulative[[1]]))
            return(as.matrix(x)[design$subsets$cumulative[[phase]],,drop=FALSE])
        stop("explicit x must match first or last phase sample size")
    }
    as.matrix(x)[design$subsets$cumulative[[phase]],,drop=FALSE]  
}


dim.multiphase<-function(x,...) dim(x$data)
dimnames.multiphase<-function(x,...) dimnames(x$data)


"[.multiphase"<-function (x, i, ..., drop = TRUE) 
{
    if (!missing(i)) {
            if (is.logical(i)) 
                x$subpop[!i] <- FALSE
            else if (is.numeric(i) && length(i)) 
                x$subpop[-i] <- FALSE
    }    else {
        if (!is.null(x$data)) 
            x$data <- x$data[, ..1, drop = FALSE]
    }
    x
}
