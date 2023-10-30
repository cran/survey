svycontrast.svyvar<-function(stat,contrasts,add=FALSE, ...){
    s<-as.vector(as.matrix(stat))
    nms<-as.vector(outer(rownames(stat),colnames(stat),paste,sep=":"))
    v<-vcov(stat)
    names(s)<-nms
    dimnames(v)<-list(nms,nms)
    attr(s,"var")<-v
    attr(s,"statistic")<-"variance"
    class(s)<-"svystat"
    svycontrast(s,contrasts=contrasts,add=add,...)
}
