\name{withCrossval}
\alias{withCrossval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Crossvalidation using replicate weights
}
\description{
  In each set of replicate weights there will be some clusters that have
  essentially zero weight. These are used as the test set, with the
  other clusters used as the training set. Jackknife weights
  ("JK1","JKn") are very similar to cross-validation at the cluster
  level; bootstrap weights are similar to bootstrapping for cross-validation.
}
\usage{
withCrossval(design, formula, trainfun, testfun, loss = c("MSE",
"entropy", "AbsError"), intercept, tuning,  nearly_zero=1e-4,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{
    A survey design object (currently only \code{svyrep.design})
  }
  \item{formula}{
    Model formula where the left-hand side specifies the outcome
    vairable and the right-hand side specifies the variables that will
    be used for prediction
  }
  \item{trainfun}{
    Function taking a predictor matrix \code{X}, an outcome vector
    \code{y},  a weights vector \code{w}, and an element of
    \code{tuning}, and training a model that is
    returned as some R object.
  }
  \item{testfun}{
    Function taking a predictor matrix \code{X} and the output from
    \code{trainfun} and returning fitted values for the outcome variable.
  }
  \item{loss}{
    Loss function for assessing prediction
  }
  \item{intercept}{
    Should the predictor matrix have an intercept added?
  }
  \item{tuning}{vector of tuning parameters, such as the regularisation
    parameter in information criteria or the number of predictors. \code{trainfun} and
    \code{testfun} will be called with each element of this vector in
    turn. Use any single-element vector if no tuning parameter is
    needed}
  \item{nearly_zero}{test-set threshold on the scale of replicate weight
  divided by sampling weight.}
  \item{\dots}{
    future expansion
  }
}

\value{
A number
}
\references{
Iparragirre, A., Lumley, T., Barrio, I., & Arostegui, I. (2023).
Variable selection with LASSO regression for complex survey data.
Stat, 12(1), e578. 
}


\seealso{
\code{\link{as.svrepdesign}}
}
\examples{
data(api)
rclus1<-as.svrepdesign(svydesign(id=~dnum, weights=~pw, data=apiclus1,
fpc=~fpc))


withCrossval(rclus1, api00~api99+ell+stype,
  trainfun=function(X,y,w,tuning) lm.wfit(X,y,w),
  testfun=function(X, trainfit,tuning) X\%*\%coef(trainfit),
  intercept=TRUE,loss="MSE",tuning=1)


## More realistic example using lasso
## tuning parameter is number of variables in model
##
##  library(glmnet)
##  ftrain=function(X,y,w,tuning) {
##   	m<-glmnet(X,y,weights=w)
##   	lambda<-m$lambda[min(which(m$df>=tuning))]
## 	list(m,lambda)
##   	}
##   ftest=function(X, trainfit, tuning){
##   	predict(trainfit[[1]], newx=X, s=trainfit[[2]])
##   }
##
##  withCrossval(rclus1, api00~api99+ell+stype+mobility+enroll,
##    trainfun=ftrain,
##    testfun=ftest,
##    intercept=FALSE,loss="MSE",
##    tuning=0:3)
##
## [1] 11445.2379  9649.1150   800.0742   787.4171
##
## Models with two or three predictors are about equally good


}