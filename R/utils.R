
#' Internal EnsBagg functions
#' 
#' @description Internal EnsBagg helper functions
#'
#' @details These functions are not intended for use by users.
#'
#' @name EnsBagg-internal

NULL

#' @rdname EnsBagg-internal
#' @param lambda penalization term. It is a positive scalar.
#' @param data A data frame that contains at least: ttilde, delta, wts
#' @param Z a matrix that contains the predictions. Each column represents a single marker.
#' @return a vector with the optimal AUC value and the optimal coefficient  


optimun_auc_coef = function(lambda,data,Z){
  optimal_coef <- optim(par = coef_init.1, fn = .cvAUC,lambda=lambda ,Z=Z, data=data,method = "BFGS",control=list(maxit=10000))
  AUC_coef_opt<- AUC_function(T=data$ttilde,delta=data$delta,marker=crossprod(t(Z),optimal_coef$par),cause=1,wts=data$wts.nn,tao)
  return(c(AUC_coef_opt,optimal_coef$convergence,optimal_coef$par))
}


#' @description  Compute the risk of missclassifying an individual using as a marker a single prediction or weighted linear combination of several predictions (1-AUC) 
#' @param par  a vector of coefficients/weights. Its length must be equal to the number of predictions included in Z 
#' @param lambda penalization term. It is a positive scalar.
#' @param Z a matrix that contains the predictions. Each column represents a single marker.
#' @param data A data frame  that constains at least: ttilde= time to event, delta=event type, wts= IPC weights
#' @return 1-AUC
#' @rdname EnsBagg-internal


risk_auc<- function(par,lambda,Z,data){
  par <- par/sum(par) 
  marker_pred <- crossprod(t(Z),par)
  Risk_AUC=1-AUC_function(T=data$ttilde,delta=data$delta,marker=marker_pred,cause=1,wts=data$wts.nn,tao)+lambda*sum(abs(par))
  return(Risk_AUC)
}

