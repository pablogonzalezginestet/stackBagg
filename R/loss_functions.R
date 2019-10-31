
#' IPC Weighted AUC Loss Function
#' @description Compute time-varying IPCW AUC to account for censoring and competing risks.
#' @param T  vector of (censored) event-times
#' @param delta  vector of event indicators at the corresponding value of the vector T. Censored observations must be denoted by the value 0. 
#' @param marker   vector of the marker values for which we want to compute the time-dependent ROC curve. the function assumes that larger values of the marker are associated with higher risks of events
#' @param cause  value of the event indicator (the non-censored observation) that represents the event of interest for which we aim to compute the time-dependent ROC curve.
#' @param wts IPC weights
#' @param tao evaluation time point of interest
#' @return value of the AUC at tao
#' @rdname ipcw_auc
#' @export

ipcw_auc <- function (T,delta, marker,cause, wts, tao){
  
  AireTrap <- function(Abs, Ord) {
    nobs <- length(Abs)
    dAbs <- Abs[-1] - Abs[-nobs]
    mil <- (Ord[-nobs] + Ord[-1])/2
    area <- sum(dAbs * mil)
    return(area)
  }
  
  n <- length(T)
  n_marker <- length(unique(marker))
  n_times <- length(tao)
  AUC <- rep(NA, n_times)
  order_marker <- order(-marker)
  Mat_data <- cbind(T,delta, marker)[order_marker, ]
  colnames(Mat_data) <- c("T","delta", "marker")
  Weights_cases_all <- wts
  Weights_cases_all <- Weights_cases_all[order_marker]
  
  Cases <- (Mat_data[, "T"] < tao & Mat_data[, "delta"] == cause)
  Controls_1 <- (Mat_data[, "T"] > tao)
  
  Weights_controls_1 <- wts
  Weights_controls_1 <- Weights_controls_1[order_marker]
  Weights_cases <- Weights_cases_all
  
  Weights_cases[!Cases] <- 0
  Weights_controls_1[!Controls_1] <- 0
  
  
  den_TP_t <- sum(Weights_cases)
  TP_tbis <- c(0, cumsum(Weights_cases))/den_TP_t
  TP_t <- TP_tbis[!duplicated(marker[order_marker])]
  
  Controls_2 <- (Mat_data[, "T"] < tao & Mat_data[,"delta"] != cause & Mat_data[, "delta"] != 0)
  Weights_controls_2 <- Weights_cases_all
  Weights_controls_2[!Controls_2] <- 0
  
  den_FP_2_t <- sum(Weights_controls_2) + sum(Weights_controls_1)
  FP_2_tbis <- c(0, cumsum(Weights_controls_1)+cumsum(Weights_controls_2) )/den_FP_2_t
  FP_2_t <- FP_2_tbis[!duplicated(marker[order_marker])]
  
  AUC <- AireTrap(FP_2_t, TP_t)
  
  return(AUC)
  
}




#' IPC Weighted Brier Loss Function
#' @description Compute weighted Brier Loss function for a single marker or a linear weighted combination of markers
#' @param par a vector of weights. Its length must be equal to the number of predictions included in Z 
#' @param Z  a matrix that contains the predictions. Each column represents a single marker.
#' @param y vector of response variable (binary).
#' @param wts IPC weights
#' @rdname wBrierScore
#' @export

ipcw_brier<- function(par,Z,y,wts){
  no.na <- !is.na(y)
  y <- as.numeric(y)
  brier <- sum((wts[no.na] * (y[no.na]-crossprod(t(Z),par)[no.na])^2))/sum(wts)
  return(brier)
}



#' IPC Weighted Cross Entropy Loss Function
#' @description Compute weighted Cross Entropy Loss function for a single marker or a linear weighted combination of markers
#' @param par a vector of weights. Its length must be equal to the number of predictions included in Z 
#' @param Z  a matrix that contains the predictions. Each column represents a single marker.
#' @param y vector of response variable (binary).
#' @param wts IPC weights
#' @rdname wCrossEntropy
#' @export


ipcw_crossentropy <- function(par,Z,y,wts){
  no.na <- !is.na(y)
  y <- as.numeric(y)-1
  eps <- 1e-15
  y_hat <- crossprod(t(Z),par)
  y_hat <- pmax(pmin(y_hat, 1 - eps), eps)
  H <- - sum(wts[no.na] * y[no.na] * log(y_hat[no.na]) + wts[no.na] * (1-y[no.na]) * log(1-y_hat[no.na]))/(sum(wts))
  return(H)
}






