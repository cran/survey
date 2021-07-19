## quantiles with equal weights

library(survey)
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)


for(i in 1:9){
    print(i)
    all.equal(
        as.vector(coef(svyquantile(~ell, dclus1, c(0.2,0.5,0.9), qrule=paste0("hf",i)))),
        as.vector(quantile(apiclus1$ell,  c(0.2,0.5,0.9), type=i))
    )
    
}
