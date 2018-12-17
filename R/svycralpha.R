
svycralpha<-function(formula, design, na.rm=FALSE){

    scoredef<-formula
    scoredef[[1]]<-quote(I)    
    design<-eval(bquote(update(design, `*alpha*`= .(scoredef))))
    vtotal<-coef(svyvar(~`*alpha*`,design,na.rm=na.rm))
    vitems <-diag(coef(svyvar(formula, design,na.rm=na.rm)))
    K<-length(attr(terms(formula),"term.labels"))
    (K/(K-1))*(1-sum(vitems)/vtotal)
}
