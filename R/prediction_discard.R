

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
                                tuneparams=NULL ) {
  
  
  xnam.factor <- colnames(train.data[xnam])[sapply(train.data[xnam], class)=="factor"]
  if(length(xnam.factor)==0){ xnam.factor<- NULL}
  xnam.cont <- xnam[!(xnam %in% xnam.factor)]
  xnam.cont.gam <- xnam.cont[apply(train.data[xnam.cont],2, function(z) length(unique(z))>3 )]
  
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
    tuneparams <- ensBagg::parametersSimulation(folds = 5,xnam,train.data,tao)
  }  
  
  if(length(tuneparams$gam)>1){
    A<- length(ML_list)+1 # if we consider gam with df=3 and df=4
    A_native <- length(ML_list_natively)+1
    ml_names<- c("LogisticReg","GAM.3","GAM.4","LASSO","Random Forest","SVM","BART","k-NN","Neural Network")
  }else{
    A<- length(ML_list)
    A_native <- length(ML_list_natively)
    ml_names<- c("LogisticReg","GAM","LASSO","Random Forest","SVM","BART","k-NN","Neural Network")
  }

  pred1 <- ML_list$logfun(train.data,test.data,fmla, xnam,xnam.factor,xnam.cont)
  pred2 <- ML_list$GAMfun(train.data,test.data,fmla,xnam,xnam.factor,xnam.cont,xnam.cont.gam,tuneparams$gam_param)
  pred3 <- ML_list$lassofun(train.data,test.data,fmla,xnam,xnam.factor,xnam.cont,tuneparams$lasso_param)
  pred4 <- ML_list$rffun(train.data,test.data,fmla,xnam,xnam.factor,xnam.cont,tuneparams$randomforest_param)
  pred5 <- ML_list$svmfun(train.data,test.data,xnam,xnam.factor,xnam.cont,tuneparams$svm_param)
  pred6 <- ML_list$bartfun(train.data,test.data,xnam,xnam.factor,xnam.cont, tuneparams$bart_param)
  pred7 <- ML_list$knnfun(train.data,test.data,xnam,xnam.factor,xnam.cont, tuneparams$knn_param)
  pred8 <- ML_list$nn(train.data,test.data,xnam,xnam.factor,xnam.cont,tuneparams$nn_param)
  
  prediction_discard<- as.matrix(cbind(pred1,pred2,pred3,pred4,pred5,pred6,pred7,pred8))
  colnames(prediction_discard) <- ml_names
  auc_discard <- apply(prediction_discard,2, function(x) cvAUC::AUC(predictions = x,labels =test.data$E) )
  
  return( list(prediction_discard=prediction_discard,auc_discard=auc_discard) )
  
}