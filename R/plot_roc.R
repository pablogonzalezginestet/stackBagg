
#' Plot IPCW ROC curve
#' @description Plot IPCW ROC curve for a prediction 
#' @param time vector of (censored) event-times
#' @param delta vector of event indicators at the corresponding value of the vector time. Censored observations must be denoted by the value 0. 
#' @param marker vector of the marker values for which we want to compute the time-dependent ROC curve. the function assumes that larger values of the marker are associated with higher risks of events
#' @param wts IPC weights
#' @param tao evaluation time point of interest 
#' @param method IPCW or discard 
#' @return a plot IPCW ROC curve
#' @import ggplot2
#' @export



plot_roc <- function(time,delta, marker, wts=NULL, tao, method){
    
    if(method=="discard"){
    E=as.factor(ifelse(time < tao & delta==1, 1 , ifelse(time < tao & delta==2 | time>tao, 0, NA) ) )
    #time <- time[!is.na(E)]
    #delta <- delta[!is.na(E)]
    marker <- marker[!is.na(E)]
    #wts <- rep(1,length(time))
    discard_roc_plot<- plot(AUC::roc(predictions=marker,labels=E[!is.na(E)]))
    return(discard_roc_plot)
    }
    
    
    if(method=="ipcw"){
    
    cause <- 1
    n <- length(time)
    n_marker <- length(unique(marker))
    n_times <- length(tao)
    AUC <- rep(NA, n_times)
    order_marker <- order(-marker)
    Mat_data <- cbind(time,delta, marker)[order_marker, ]
    colnames(Mat_data) <- c("time","delta", "marker")
    Weights_cases_all <- wts
    Weights_cases_all <- Weights_cases_all[order_marker]
    
    Cases <- (Mat_data[, "time"] < tao & Mat_data[, "delta"] == cause)
    Controls_1 <- (Mat_data[, "time"] > tao)
    
    Weights_controls_1 <- wts
    Weights_controls_1 <- Weights_controls_1[order_marker]
    Weights_cases <- Weights_cases_all
    
    Weights_cases[!Cases] <- 0
    Weights_controls_1[!Controls_1] <- 0
    
    
    den_TP_t <- sum(Weights_cases)
    TP_tbis <- c(0, cumsum(Weights_cases))/den_TP_t
    TP_t <- TP_tbis[!duplicated(marker[order_marker])]
    
    Controls_2 <- (Mat_data[, "time"] < tao & Mat_data[,"delta"] != cause & Mat_data[, "delta"] != 0)
    Weights_controls_2 <- Weights_cases_all
    Weights_controls_2[!Controls_2] <- 0
    
    den_FP_2_t <- sum(Weights_controls_2) + sum(Weights_controls_1)
    FP_2_tbis <- c(0, cumsum(Weights_controls_1)+cumsum(Weights_controls_2) )/den_FP_2_t
    FP_2_t <- FP_2_tbis[!duplicated(marker[order_marker])]
    
    auc_temp <- as.data.frame(cbind(y=TP_t,x=FP_2_t))
    
    # This is the plot!!!
    ipcw_roc_plot <- ggplot(auc_temp, aes(x = x, y = y)) +
    geom_line()+
    xlab("\n False Negative (1 - Sp)")+
    ylab("True Positive (Se) \n")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),plot.margin = margin(10,0,10,0),
    axis.title.x = element_text( size=16),
    axis.title.y = element_text( size=16),
    axis.text.x = element_text( size=14),
    axis.text.y = element_text( size=14))+
    theme(aspect.ratio = 1)+
    geom_abline(slope =1,linetype="dotted")
    
    
    return(ipcw_roc_plot)
    
    }
    
    }

