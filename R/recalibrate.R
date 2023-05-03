
recalibrate<-function(design, formula, ...){

    m<-model.matrix(formula, model.frame(design))
    calibrate(design,formula, colSums(m*weights(design)))
    
}
