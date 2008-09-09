svysmooth<-function(formula,design,method="locpoly",bandwidth,...) UseMethod("svysmooth", design)
svysmooth.default<-function(formula, design,method="locpoly",bandwidth,...){
  switch(match.arg(method),
         locpoly=svylocpoly(formula,design,bandwidth=bandwidth,...))
}

svylocpoly<-function(formula, design, ngrid=401, xlim=NULL,
                     ylim=NULL, bandwidth,...){
  require("KernSmooth")

  mf<-model.frame(formula,model.frame(design))
  mm<-model.matrix(terms(formula),mf)
  if(attr(terms(formula),"intercept"))
    mm<-mm[,-1,drop=FALSE]

  bandwidth<-rep(bandwidth, length=ncol(mm))

  if (length(formula)==3){
    Y<-model.response(mf)
    density<-FALSE
  } else density<-TRUE
  

  if (is.null(xlim)){
    xlim<-apply(mm,2,range)
  }
  if (!is.matrix(xlim))
    xlim<-matrix(xlim,nrow=2)
  
  w<-weights(design,type="sampling")

  ll<-vector("list", ncol(mm))
  for(i in 1:NCOL(mm)){
    gx<-seq(min(xlim[,i]), max(xlim[,i]), length=ngrid)
    nx<-rowsum(c(rep(0,ngrid),w), c(1:ngrid, findInterval(mm[,i],gx)))
    if (density){
      ll[[i]]<-locpoly(rep(1,ngrid),nx*ngrid/(diff(xlim[,i])*sum(w)),
                           binned=TRUE, bandwidth=bandwidth[i], range.x=xlim[,i])
    }else{
      ny<-rowsum(c(rep(0,ngrid), Y*w), c(1:ngrid, findInterval(mm[,i],gx)))
      ll[[i]]<-locpoly(nx, ny, binned=TRUE, bandwidth=bandwidth[i], range.x=xlim[,i])
    }
    names(ll)<-attr(terms(formula),"term.labels")
  }
  attr(ll,"call")<-sys.call(-2)
  attr(ll,"density")<-density
  if(density)
    attr(ll,"ylab")<-"Density"
  else
    attr(ll,"ylab")<-deparse(formula[[2]])

  class(ll)<-"svysmooth"
  
  ll
 
}

print.svysmooth<-function(x,...){
  if(attr(x,"density"))
    cat("Density estimate: :")
  else
    cat("Scatterplot smoother :")
  print(attr(x,"call"))
  invisible(x)
}

plot.svysmooth<-function(x, which=NULL,type="l",xlabs=NULL,ylab=NULL,...){
  if (is.null(which))
    which<-seq(length=length(x))
  if (is.character(which))
    which<-match(which,names(x))
  
  if(is.null(xlabs)) xlabs<-names(x)[which]
  if(is.null(ylab)) ylab<-attr(x,"ylab")

  for(i in seq(length=length(which)))
    plot(x[[which[i]]], type=type, xlab=xlabs[i], ylab=ylab, ...)
  invisible(NULL)
}

lines.svysmooth<-function(x,which=NULL,...){
  for(i in names(x)) lines(x[[i]],...)
}
