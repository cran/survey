##
## tests calibration using a working model 
##
library(survey)
data(nwtco)
 nwtco$incc2<-as.logical(with(nwtco, ifelse(rel | instit==2,1,rbinom(nrow(nwtco),1,.1))))

dccs2<-twophase(id=list(~seqno,~seqno),strata=list(NULL,~interaction(rel,instit)),
    data=nwtco, subset=~incc2)
    
impmodel<-svyglm(I(histol==2)~instit*stage+age, design=dccs2, family=quasibinomial())
nwtco$histhat<-1+predict(impmodel,newdata=nwtco)

phase1model<-glm(rel~histhat*factor(stage),data=nwtco)

newmethod<- calibrate(dccs2, formula=phase1model,phase=2,calfun=cal.raking)

infs<-survey:::estfuns(phase1model)
colnames(infs)<-paste0("h",1:ncol(infs))
nwtco<-cbind(nwtco,infs)

dccs2h<- twophase(id=list(~seqno,~seqno),strata=list(NULL,~interaction(rel,instit)),
    data=nwtco, subset=~incc2)
oldmethod<-calibrate(dccs2h, ~h1+h2+h3+h4+h5+h6+h7+h8,calfun="raking",phase=2)

a<-svyglm(I(histol==2)~instit*stage+age, design=newmethod, family=quasibinomial())
b<-svyglm(I(histol==2)~instit*stage+age, design=oldmethod, family=quasibinomial())

all.equal(coef(a),coef(b))
all.equal(SE(a), SE(b),tol=1e-3)
