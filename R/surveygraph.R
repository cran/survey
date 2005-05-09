
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
           } else {
             ## new version
             rval<-hexbin(X,Y,ID=TRUE)
             cell<-rval@cID
             rval@count<-as.vector(tapply(W,cell,sum))
             rval@xcm<-as.vector(tapply(1:length(X), cell,
                              function(ii) weighted.mean(X[ii],W[ii])))
             rval@ycm<-as.vector(tapply(1:length(Y), cell,
                              function(ii) weighted.mean(Y[ii],W[ii])))
           }
           
           plot(rval,legend=legend,style="centroids",...)
         },
         grayhex={
           if (exists("hcell")) {
             ## old version of hexbin
             rval<-hexbin(X,Y)
             cell<-hcell(X,Y)$cell
             rval$cnt<-tapply(W,cell,sum)
         } else {
             ## new version
             rval<-hexbin(X,Y,ID=TRUE)
             cell<-rval@cID
             rval@count<-as.vector(tapply(W,cell,sum))
           }
           plot(rval, legend=legend,...)
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
