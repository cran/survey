#' calculations for Bell-McCaffrey standard errors and degrees of freedom
#' 
#' @param resid residuals (vector) or estimating functions (matrix)
#' @param design the survey design, trimmed to missing data and subset conditions
#' @param regressors design matrix of the glm predictors from model.matrix()
#' @param std.errors "Bell-McCaffrey" for Bell-McCaffrey standard errors 
#' with unit working covariance matrix; "
#' Bell-McCaffrey-2" for Bell-McCaffrey standard errors 
#' with exchangeable correlation working covariance matrix
#' @param degf whether the degrees of freedom are to be computed
#' @returns the vector of scaled residuals
#
# does not have to be exported
scale_bell_mcaffrey<-function(resid,design,regressors,std.errors, degf=FALSE) {
  
  if (!is.null(getOption("svy.debug.bmca"))) browser()
  
  resid2<-as.matrix(resid)
  
  if (nrow(resid2)!=nrow(design$variables) | nrow(resid2)!=nrow(regressors))
    stop("length mismatch: this should not be happening.")
  
  XtX<-crossprod(regressors)
  invXtX<-chol2inv(chol(XtX))
  
  if (degf) {
    # nrow(data) x nrow(data) matrix -- should it be sparse??
    H <- tcrossprod(tcrossprod(regressors, invXtX), regressors)
    H <- diag(nrow(design$variables))-H
    G <- list()
  }
  
  for (k in 1:length(unique(design$cluster[,1]))) {
    # pick up the cluster ID
    clID<-unique(design$cluster)[k,1]
    # pick up the row indices
    cl_indices<-which(design$cluster[,1]==clID)
    
    # if (length(cl_indices)>ncol(regressors)) {
    # slice the regressor matrix; do not let R convert it to a vector!
    this_X<-regressors[cl_indices,,drop=FALSE]
    # slice the residual matrix
    this_res<-resid2[cl_indices,]
    # form the projection
    this_proj<-tryCatch( tcrossprod( tcrossprod(this_X, invXtX), this_X), 
                         error = function(e) {
                           cat("Loop index = ", k, "\n")
                           cat("Cluster ID = ",clID, "\n")
                           cat("Cluster indices = ", cl_indices, "\n")
                           cat("Regressor matrix dimensions are: \n", dim(this_X), "\n")
                           cat("(X'X)^{-1} matrix dimenstions are:\n", dim(invXtX), "\n")
                           cat("Residual vector is:\n", this_res, "\n")
                           # return zero matrix if the projection fails
                           matrix(rep(0,length(cl_indices)*length(cl_indices)), nrow=length(cl_indices))
                         }
                         
    )
    # update the residuals
    A_i <- mihalf(diag(length(cl_indices))-this_proj)
    resid2[cl_indices,]<-crossprod(A_i, this_res)
    
    # for the degrees of freedom calculation
    if (degf) {
      # Theorem 4 of Bell-McCaffrey (SMJ 2002);
      # this is the transpose of g_i
      # the dimension is {nrow(design$variables)} x {# of regressors} 
      # dfs are specific to parameters and have to be extracted by component
      G[[k]] <- tcrossprod(crossprod(H[cl_indices,,drop=FALSE], A_i), tcrossprod( invXtX, this_X ))
    }
    # } endif dim regressors
  }
  
  if (degf) {
    if (std.errors=="Bell-McCaffrey-2") {
      # exchangeable correlation matrix, form an estimate
      sigma2 <- sum(resid2*resid2)/length(resid2)
      sigma_cross <- 0
      n_cross <- 0
      for (k in 1:length(unique(design$cluster[,1]))) {
        # pick up the cluster ID
        clID<-unique(design$cluster)[k,1]
        # pick up the row indices
        cl_indices<-which(design$cluster[,1]==clID)
        
        sigma_cross <- sigma_cross + sum(tcrossprod(resid2[cl_indices]))
        n_cross <- n_cross + length(cl_indices)*length(cl_indices)
      }
      sigma_cross <- (sigma_cross - sum(resid2*resid2)) / (n_cross - length(resid2))
      
      # form monstrous V; should be a sparse matrix
      V <- diag(sigma2,nrow=nrow(design$variables))
      # replace off-diagonal elements within the same cluster
      for (k in 1:length(unique(design$cluster[,1]))) {
        # pick up the cluster ID
        clID<-unique(design$cluster)[k,1]
        # pick up the row indices
        cl_indices<-which(design$cluster[,1]==clID)
        V[cl_indices, cl_indices] <- sigma_cross
      }
      # the diagonal entries got overwritten
      diag(V) <- rep(sigma2, nrow(V))
    }
    
    lambda <- numeric()
    for (j in 1:ncol(regressors)) {
      # extract the j-th component from each G[[k]] 
      this_G <- sapply(G, function(X, col) X[,col], j)
      # this is supposed to be nrow(design$variables) x {# clusters} matrix
      if (std.errors=="Bell-McCaffrey") {
        # identify matrix, simple!!!
        GVG <- tcrossprod(this_G)
      } else if (std.errors=="Bell-McCaffrey-2") {
        # exchangeable correlation matrix
        GVG <- crossprod(crossprod(V, this_G), this_G)
      } else {
        stop("Standard errors type not supported: ", std.errors)
      }
      eigen_GVG <- eigen(GVG, symmetric=TRUE, only.values=TRUE)
      this_df <- { sum(eigen_GVG$values)*sum(eigen_GVG$values) / 
          sum(eigen_GVG$values*eigen_GVG$values) }
      lambda[j]<-this_df
    }
    attr(resid2,"degf") <- lambda
  }
  
  return(resid2)
}

#' Symmetric square root of a matrix
#' 
#' Computes the symmetric square root of a positive definite matrix
#' 
#' 
#' @usage mhalf(M)
#' @param M a positive definite matrix
#' @return a matrix \code{H} such that \code{H^2} equals \code{M}
#' @author Peter Hoff
#' @export 
mhalf <- function(M) { 
    #symmetric square  root of a pos def matrix
    tmp<-eigen(M)
    tmp$vec%*%sqrt(diag(tmp$val,nrow=nrow(M)))%*%t(tmp$vec)
  }

#' Symmetric inverse square root of a matrix
#' 
#' Computes the inverse symmetric square root of a positive definite matrix
#' 
#' 
#' @usage mihalf(M)
#' @param M a positive definite matrix
#' @return a matrix $H$ such that $H^2$ equals $M^{-1]$
#' @author Peter Hoff, Stas Kolenikov
#' @export 
mihalf <- function(M) { 
    #symmetric inverse square root of a pos def matrix
    tmp<-eigen(M)
    eigen0<-1/tmp$val
    eigen0[tmp$val==0]<-0
    tmp$vec%*%sqrt(diag(eigen0,nrow=nrow(M)))%*%t(tmp$vec)
  }
