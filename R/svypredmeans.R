
svypredmeans<-function(adjustmodel, groupfactor){

	design<-eval(bquote(update(adjustmodel$survey.design, .groupfactor=.(groupfactor[[2]]))))
	groups<-unique(model.frame(design)$.groupfactor)
	groups<-groups[!is.na(groups)]
	model<-update(adjustmodel, .~.+.groupfactor,design=design)
	w<-weights(design,"sampling")
	
	fits<-matrix(nrow=NROW(design),ncol=length(groups))
	dg_deta<-matrix(nrow=length(coef(model)),ncol=length(groups))
	for(i in 1:length(groups)){
		mf<-model.frame(design)
		mf$.groupfactor<-groups[i]
		mu<-predict(model,newdata=mf,type="response",se.fit=FALSE)
		eta<-predict(model,newdata=mf,type="link",se.fit=FALSE)
		fits[,i]<-coef(mu)
		
		mm<-model.matrix(terms(model),mf)
		dg_deta[,i]<-t(colSums(w*model$family$mu.eta(eta)*mm))/sum(w)
	}
	colnames(fits)<-as.character(groups)
	cond<-svymean(fits,design)
	addvar<-t(dg_deta)%*%vcov(model)%*%dg_deta
	vv<-addvar+attr(cond,"var")
	attr(vv,"parts")<-list(addvar,attr(cond,"var"))
	attr(cond,"var")<-vv
	cond
}

