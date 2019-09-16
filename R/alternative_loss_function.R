
#' Weighted Brier Loss Function
#' @description Compute weighted Brier Loss function for a single marker or a linear weighted combination of markers
#' @param par a vector of weights. Its length must be equal to the number of predictions included in Z 
#' @param Z  a matrix that contains the predictions. Each column represents a single marker.
#' @param y vector of response variable (binary).
#' @param wts IPC weights
#' @rdname wBrierScore
#' @export

wBrierScore<- function(par,Z,y,wts){
  no.na <- !is.na(y)
  y <- as.numeric(y)
  brier <- sum((wts[no.na] * (y[no.na]-crossprod(t(Z),par)[no.na])^2))/sum(wts)
  return(brier)
}



#' Weighted Cross Entropy Loss Function
#' @description Compute weighted Cross Entropy Loss function for a single marker or a linear weighted combination of markers
#' @param par a vector of weights. Its length must be equal to the number of predictions included in Z 
#' @param Z  a matrix that contains the predictions. Each column represents a single marker.
#' @param y vector of response variable (binary).
#' @param wts IPC weights
#' @rdname wCrossEntropy
#' @export


wCrossEntropy <- function(par,Z,y,wts){
  no.na <- !is.na(y)
  y <- as.numeric(y)-1
  eps <- 1e-15
  y_hat <- crossprod(t(Z),par)
  y_hat <- pmax(pmin(y_hat, 1 - eps), eps)
  H <- - sum(wts[no.na] * y[no.na] * log(y_hat[no.na]) + wts[no.na] * (1-y[no.na]) * log(1-y_hat[no.na]))/(sum(wts))
  return(H)
}






