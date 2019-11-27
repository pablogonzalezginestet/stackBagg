

#' Prediction discarding censored observations
#' @description Ad-hoc technique that discards censored observations from the analysis 
#' @param train.data a data.frame with at least the following variables: 
#' event-times (censored) in the first column,
#' event indicator in the second column and 
#' covariates/features that the user potentially want to use in building the preodiction model.
#' Censored observations must be denoted by the value 0. Main event of interest is denoted by 1.
#' @param test.data a data.frame with the same variables and names that the train.data  
#' @param xnam vector with the names of the covariates to be included in the model
#' @param tao evaluation time point of interest
#' @param tuneparams a list of tune parameters for each machine learning procedure. Name them as gam_param, lasso_param, randomforest_param, svm_param, bart_param, knn_param, nn_param.
#' Default values are the same used for the simulation.
#' @return a list with the predictions of each machine learning algorithm and the AUC of each of them
#' @export

prediction_discard <- function( train.data,
                                test.data ,
                                xnam,
                                tao,
                                ens.library,
                                tuneparams=NULL ) {
  
  
  xnam.factor <- colnames(train.data[xnam])[sapply(train.data[xnam], class)=="factor"]
  if(length(xnam.factor)==0){ xnam.factor<- NULL}
  xnam.cont <- xnam[!(xnam %in% xnam.factor)]
  xnam.cont.gam <- xnam.cont[apply(train.data[xnam.cont],2, function(z) length(unique(z))>3 )]
  
  all.library <- stackBagg::all.algorithms()
  
  names(train.data)[1] <- "ttilde"
  names(train.data)[2] <- "delta"
  #names(test.data)[1] <- "id"
  names(test.data)[1] <- "ttilde"
  names(test.data)[2] <- "delta"
  
  train.data<- dplyr::mutate(train.data,id = 1:length(ttilde),E=as.factor(ifelse(ttilde < tao & delta==1, 1 , ifelse(ttilde < tao & delta==2 | ttilde>tao, 0, NA))),
                             deltac=ifelse(delta==0,1,0))
  test.data<- dplyr::mutate(test.data,id = 1:length(ttilde),E=as.factor(ifelse(ttilde < tao & delta==1, 1 , ifelse(ttilde < tao & delta==2 | ttilde>tao, 0, NA))),
                            deltac=ifelse(delta==0,1,0))
  
  train.data$wts <- 1
  test.data$wts <- 1
  train.data <- train.data[c("id","E","wts","ttilde","delta",xnam)]
  test.data <- test.data[c("id","E","wts","ttilde","delta",xnam)]
  
  train.data <- train.data[!is.na(train.data$E),]
  test.data <- test.data[!is.na(test.data$E),]
  
  
  fmla <- as.formula(paste("E ~ ", paste(xnam, collapse= "+")))
  
  if (missing(tuneparams)) {
    tuneparams <- stackBagg::parametersSimulation(folds = 5,xnam,train.data,tao)
  }  
  
  
  if("ens.gam" %in% ens.library & length(tuneparams$gam)>1){
    ml_names<- all.library[all.library %in% ens.library]
    ml_names <- c(ml_names[1],"ens.gam.3","ens.gam.4",ml_names[-(1:2)])
   
  }else{
    ml_names<- all.library[all.library %in% ens.library]
  }
  
  pred <- NULL
  if("ens.glm" %in% ens.library){
  pred_temp <- ML_list$logfun(train.data,test.data,fmla, xnam,xnam.factor,xnam.cont)
  pred <- cbind(pred_temp,pred)
  }
  if("ens.gam" %in% ens.library){
  pred_temp <- ML_list$GAMfun(train.data,test.data,fmla,xnam,xnam.factor,xnam.cont,xnam.cont.gam,tuneparams$gam_param)
  pred <- cbind(pred_temp,pred)
   }
  if("ens.lasso" %in% ens.library){
  pred_temp <- ML_list$lassofun(train.data,test.data,fmla,xnam,xnam.factor,xnam.cont,tuneparams$lasso_param)
  pred <- cbind(pred_temp,pred)
  }
  if("ens.randomForest" %in% ens.library){
  pred_temp <- ML_list$rffun(train.data,test.data,fmla,xnam,xnam.factor,xnam.cont,tuneparams$randomforest_param)
  pred <- cbind(pred_temp,pred)
  }
  if("ens.svm" %in% ens.library){
  pred_temp <- ML_list$svmfun(train.data,test.data,xnam,xnam.factor,xnam.cont,tuneparams$svm_param)
  pred <- cbind(pred_temp,pred)
  }
  if("ens.bartMachine" %in% ens.library){
  pred_temp <- ML_list$bartfun(train.data,test.data,xnam,xnam.factor,xnam.cont, tuneparams$bart_param)
  pred <- cbind(pred_temp,pred)
  }
  if("ens.knn" %in% ens.library){
  pred_temp <- ML_list$knnfun(train.data,test.data,xnam,xnam.factor,xnam.cont, tuneparams$knn_param)
  pred <- cbind(pred_temp,pred)
  }
  if("ens.nn" %in% ens.library){
  pred_temp <- ML_list$nn(train.data,test.data,xnam,xnam.factor,xnam.cont,tuneparams$nn_param)
  pred <- cbind(pred_temp,pred)
  }
  
  prediction_discard<- as.matrix(pred)
  colnames(prediction_discard) <- ml_names
  auc_discard <- apply(prediction_discard,2, function(x) cvAUC::AUC(predictions = x,labels =test.data$E) )
  
  return( list(prediction_discard=prediction_discard,auc_discard=auc_discard) )
  
}