
## Need to call qs() from qrule(), because definition of p varies by rule. <sigh>

## Discrete

qrule_math <-qrule_hf1 <- function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
        }
    qdata<-qs(x,w,p)
    if(qdata$wlow==0) qdata$qlow else qdata$qup
}

qrule_school<-qrule_hf2 <- function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
        }
    qdata<-qs(x,w,p)
    if(qdata$wlow==0) (qdata$qlow+qdata$qup)/2 else qdata$qup
}

qrule_hf3 <- function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    w<-rowsum(w,x)
    x<-sort(unique(x))
    qdata<-qs(x,w,p)
    if((qdata$wlow==0)&& (qdata$ilow %%2 ==0)) qdata$qlow else qdata$qup
}

## Continuous

qrule_hf4 <- function(x,w,p){
     if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
     }
     qdata<-qs(x,w,p)
     gamma<-with(qdata, wlow/(wup+wlow))
     qdata$qlow*(1-gamma)+qdata$qup*gamma
}



qrule_hf5<-function(x,w,p){
     if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    n<-length(x)
     if (n==1) return(x)
     
     ii<-order(x)
     x<-x[ii]
     w<-w[ii]
     cumw<-cumsum(w)
     pk<-(cumw-w/2)/(cumw[n])
     approx( pk,x, p, method="linear", rule=2)$y
}



qrule_hf6<-function(x,w,p){
     if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    n<-length(x)
     if (n==1) return(x)

     ii<-order(x)
     x<-x[ii]
     w<-w[ii]
     cumw<-cumsum(w)
     pk<-cumw/(cumw[n]+w[n])
     approx( pk,x, p, method="linear", rule=2)$y
}

qrule_shahvaish<-function(x,w,p){
     if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    n<-length(x)
     if (n==1) return(x)

     ii<-order(x)
     x<-x[ii]
     w<-w[ii]
     wbar<-w/mean(w)
     S<-cumsum(wbar)
     pk<-(S+0.5-w/2)/(n+1)
     approx( pk,x, p, method="constant", f=0, rule=2)$y
}




qrule_hf7<-function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    n<-length(x)
     if (n==1) return(x)

    ii<-order(x)
    x<-x[ii]
    cumw<-cumsum(w[ii])
    pk<-c(0, cumw[-n])/cumw[n-1]
    approx( pk,x, p, method="linear", rule=2)$y
    
 }


qrule_hf8<-function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    n<-length(x)
     if (n==1) return(x)

    ii<-order(x)
    x<-x[ii]
    w<-w[ii]
    cumw<-cumsum(w)
    pk<-(c(0, cumw[-n])*1/3+cumw*2/3)/(cumw[n]+w[n]/3)
    approx( pk,x, p, method="linear", rule=2)$y
    
}


qrule_hf9<-function(x,w,p){
    if (any(zero<-w==0)){
        w<-w[!zero]
        x<-x[!zero]
    }
    n<-length(x)
     if (n==1) return(x)

    ii<-order(x)
    x<-x[ii]
    w<-w[ii]
    cumw<-cumsum(w)
    pk<-(c(0, cumw[-n])*3/8+cumw*5/8)/(cumw[n]+w[n]/4)
    approx( pk,x, p, method="linear", rule=2)$y
    
 }

last<-function(a) {
    if (any(a)) max(which(a)) else 1
    }

qs <- function(x, w, p){
    ## already has missings removed, ties handled.
    n<-length(x)
    ii<-order(x)
    x<-x[ii]
    cumw<-cumsum(w[ii])

    pos<-last(cumw<=p*sum(w))
    posnext<-if(pos==length(x)) pos else pos+1
    
    list(qlow=x[pos], qup=x[posnext], ilow=pos,iup=posnext, wlow=p-cumw[pos]/sum(w),wup=cumw[posnext]/sum(w)-p)
    
}
