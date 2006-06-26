
svyplot<-function(formula,
                  design,
                  style=c("bubble","hex","grayhex","subsample"),
                  sample.size=500, subset=NULL,legend=1,inches=0.05,...){
  
  style<-match.arg(style)
  if (style %in% c("hex","grayhex") && !require(hexbin)){
    stop(style," plots require the hexbin package")
  }

  subset<-substitute(subset)
  subset<-with(design$variables, subset)
  if(length(subset)>0)
    design<-design[subset,]

  W<-weights(design, "sampling")

  mf<-model.frame(formula, design$variables,na.action=na.pass)  
  Y<-model.response(mf)
  X<-mf[,attr(attr(mf,"terms"),"term.labels")]
  
  switch(style, 
         bubble={
           symbols(X,Y,circles=sqrt(W),inches=inches,...)
         },
         hex={
           if (exists("hcell")) {
             ## old version of hexbin
             rval<-hexbin(X,Y)
             cell<-hcell(X,Y)$cell
             rval$cnt<-tapply(W,cell,sum)
           rval$xcm<-tapply(1:length(X), cell,
                          function(ii) weighted.mean(X[ii],W[ii]))
           rval$ycm<-tapply(1:length(Y), cell,
                           function(ii) weighted.mean(Y[ii],W[ii]))
             plot(rval,legend=legend,style="centroids",...)
           } else {
             ## new version
             rval<-hexbin(X,Y,ID=TRUE)
             cell<-rval@cID
             rval@count<-as.vector(tapply(W,cell,sum))
             rval@xcm<-as.vector(tapply(1:length(X), cell,
                              function(ii) weighted.mean(X[ii],W[ii])))
             rval@ycm<-as.vector(tapply(1:length(Y), cell,
                              function(ii) weighted.mean(Y[ii],W[ii])))
             gplot.hexbin(rval, legend=legend, style="centroids",...)
           }
           
         },
         grayhex={
           if (exists("hcell")) {
             ## old version of hexbin
             rval<-hexbin(X,Y)
             cell<-hcell(X,Y)$cell
             rval$cnt<-tapply(W,cell,sum)
             plot(rval, legend=legend,...)
           } else {
             ## new version
             rval<-hexbin(X,Y,ID=TRUE)
             cell<-rval@cID
             rval@count<-as.vector(tapply(W,cell,sum))
             gplot.hexbin(rval, legend=legend,...)
           }

         },
         subsample={
           index<-sample(length(X),sample.size,replace=TRUE, prob=W)
           if (is.numeric(X))
             xs<-jitter(X[index],factor=3)
           else
             xs<-X[index]
           if (is.numeric(Y))
             ys<-jitter(Y[index],factor=3)
           else
             ys<-Y[index]
           plot(xs,ys,...)
         })

}

svyboxplot<-function(formula, design, ...){
    
    formula<-as.formula(formula)
    if(length(formula)!=3) stop("need a two-sided formula")
    if(length(formula[[3]])!=1) stop("only one rhs variable allowed")
    outcome<-eval(bquote(~.(formula[[2]])))
    
    if (length(attr(terms(formula),"term.labels"))){
        groups<-eval(bquote(~.(formula[[3]])))
        qs <- svyby(outcome,groups,design,svyquantile,ci=FALSE,
                    keep.var=FALSE,
                    quantiles=c(0,0.25,0.5,0.75,1))
        n<-NCOL(qs)
        iqr<- qs[,n-1]-qs[,n-3]
        low<-pmax(qs[,n-4],qs[n-2]-1.5*iqr)
        hi<-pmin(qs[,n],qs[n-1]+1.5*iqr)
        stats<-t(as.matrix(cbind(low,qs[,n-(3:1)],hi)))
        z<-list(stats=stats,n=coef(svytotal(groups,design)))
        for(i in 1:ncol(stats)){
            out<-c(if(qs[i,n]!=hi[i]) qs[i,n],
                   if(qs[i,n-4]!=low[i])qs[i,n-4])
            z$out<-c(z$out,out)
            z$group<-c(z$group,rep(i,length(out)))
            z$names<-as.character(qs[,1])
        }
    } else {
        qs<-svyquantile(outcome,design,ci=FALSE,
                        quantiles=c(0,0.25,0.5,0.75,1))
        iqr<-qs[4]-qs[2]
        z<-list(stats=matrix(c(max(qs[1],qs[2]-1.5*iqr),
                qs[2:4],min(qs[5],qs[4]+1.5*iqr))),
                n=sum(weights(design,"sampling")))
        z$out<-c(if(qs[5]!=z$stats[5]) qs[5],
                 if(qs[1]!=z$stats[1]) qs[1])
    }
    bxp(z,...)
}



