library(survey)
set.seed(42)

df<-data.frame(x=exp(rnorm(1000)))
df$y<-round(df$x,1)
ddf<-svydesign(id=~1,data=df)
rdf<-as.svrepdesign(ddf)

SE(oldsvyquantile(~x,ddf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE))

SE(oldsvyquantile(~x,rdf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE))

SE(oldsvyquantile(~x,ddf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE,df=Inf))

SE(oldsvyquantile(~x,rdf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE,df=Inf))


oldsvyquantile(~y,ddf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE,ties="rounded",interval.type="betaWald")

oldsvyquantile(~y,rdf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE)

oldsvyquantile(~y,ddf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE,ties="rounded",interval.type="betaWald",df=Inf)

oldsvyquantile(~y,rdf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE, df=Inf)



df<-data.frame(x=exp(rnorm(20)))
df$y<-round(df$x,1)
ddf<-svydesign(id=~1,data=df)
rdf<-as.svrepdesign(ddf)
SE(oldsvyquantile(~x,ddf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE))

SE(oldsvyquantile(~x,rdf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE))

SE(oldsvyquantile(~x,ddf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE,df=Inf))

SE(oldsvyquantile(~x,rdf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE,df=Inf))


oldsvyquantile(~y,ddf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE,ties="rounded",interval.type="betaWald")

oldsvyquantile(~y,rdf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE)

oldsvyquantile(~y,ddf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE,ties="rounded",interval.type="betaWald",df=Inf)

oldsvyquantile(~y,rdf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE, df=Inf)
