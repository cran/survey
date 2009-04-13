
marginpred<-function(model, adjustfor, predict.at, ...) UseMethod("marginpred", model)

##
## Basic strategy: calibrate on ~model*adjustfor, to set interactions to zero
##

marginpred.svycoxph<-function(model, adjustfor, predict.at,se=FALSE,...){

  if(NROW(predict.at)==0) return(NULL)
  
  design<-model$survey.design
  if (!is.null(model$na.action)) design<-design[-model$na.action,]

  modelformula<-formula(model)
  calformula<-eval(bquote( ~(.(modelformula[[3]]))*(.(adjustfor))))

  adjmf<-model.frame(terms(adjustfor), model.frame(design))
  adjmm<-model.matrix(terms(adjustfor), adjmf)
  modelmm<-model.matrix(model)[,-1,drop=FALSE]

  adjpop<-drop(colSums(adjmm*weights(design)))
  modelmeans<-drop(colSums(modelmm*weights(design))/sum(weights(design)))
  pop<-as.vector(outer(adjpop,model))
  
  
  
  
}
