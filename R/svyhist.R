svyhist<-function(formula, design, breaks = "Sturges", 
                  include.lowest = TRUE, right = TRUE, xlab=NULL,
                  main=NULL, probability=TRUE,
                  freq=!probability,...){
    mf<-model.frame(formula,design$variables, na.action=na.pass)
    if (ncol(mf)>1) stop("Only one variable allowed.")
    variable<-mf[,1]
    varname<-names(mf)
    h <- hist(variable,  plot=FALSE)
    props <- coef(svymean(~cut(variable, h$breaks),
                          design, na.rm=TRUE))
    h$density<-props/diff(h$breaks)
    h$counts <- props*sum(weights(design))
    if (is.null(xlab)) xlab<-varname
    if (is.null(main)) main<-paste("Histogram of",varname)
    plot(h, ..., freq=freq,xlab=xlab,main=main)
}
