
withPV.survey.design<-withPV.svyrep.design<-function(mapping, data, action, rewrite=TRUE,...){
	
    if(inherits(mapping,"formula")) mapping<-list(mapping)
    
    if (!is.list(mapping)) 
        stop("'mapping' must be a list of formulas")
    
    if (!all(sapply(mapping, length)==3)) 
        stop("'mapping' must be a list of two-sided formulas")
    
    df<-model.frame(data,na.action=na.pass)
    PVframes<-lapply(mapping, 
                     function(f) model.frame(f[-2], model.frame(data,na.action=na.pass)))
    nvars<-length(PVframes)
    PVnames<-sapply(mapping, function(f) deparse(f[[2]]))
    if (any(PVnames %in% colnames(data)))
        stop("working PV names must not already occur in the data")
    
    nreps<-sapply(PVframes, NCOL)
    if (length(unique(nreps))>1) 
        stop("number of plausible values must be the same for all variables")		
    nreps<-nreps[1]
    
    results<-vector("list",nreps)
    
    if(rewrite){
        sublist<-vector("list",nvars)
        names(sublist)<-PVnames
        for(i in 1:nreps){
            for(j in 1:nvars) sublist[[j]]<-as.name(names(PVframes[[j]])[i])
            
            if (is.function(action)){
                actioni<-action
                body(actioni) <- eval(bquote(substitute(.(body(actioni)), sublist)))
                results[[i]]<- actioni(data)
            } else {
                actioni <- eval(bquote(substitute(.(action), sublist)))
                results[[i]] <- eval(actioni)	
            }
        }
        
        
    } else {
        .DESIGN<-data
        
	for(i in 1:nreps){
            dfi<-lapply(PVframes, function(d) d[[i]])
		names(dfi)<-PVnames
            .DESIGN$variables<-cbind(df, as.data.frame(dfi))
            if (is.function(action))
                results[[i]] <- action(.DESIGN)
            else
                results[[i]] <- eval(action)	
	}

        }
    attr(results,"call")<-sys.call()
    results
}
    
