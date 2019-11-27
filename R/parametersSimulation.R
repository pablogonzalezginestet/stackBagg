
#' Parameter Values for the Simulation
#' @description  Values for the hyperparameters used in the simulation setup 
#' @param folds number of folds
#' @param fmla formula object ex. "E ~ x1+x2"
#' @param xnam a vector with the covariates names considered in the modeling
#' @param data a training data set
#' @return a list of values for each hyperparameter in each algorithm used in the simulation setup 
#' @export


parametersSimulation <- function(folds,xnam,data,tao){
  
  fmla <- as.formula(paste("E ~ ", paste(xnam, collapse= "+")))
  
  gam_param <- c(3,4)
  
  lasso_param <- tune_lasso(folds, fmla , tao , data, xnam )
  
  num.trees <- 500
  mtry <- round(sqrt(length(xnam))) #  the (rounded down) square root of the number variables 
  randomforest_param <- c(num.trees,mtry)
  
  knn_param <- 25
  
  cost <- .01
  gamma <- NA
  svm_param <- c(cost,gamma,2)
  
  nn_param <- 1
  
  bart_param <- c(50,1,.9)
  
  return( list(
    gam_param=gam_param,
    lasso_param=lasso_param,
    randomforest_param=randomforest_param,
    knn_param=knn_param,
    svm_param=svm_param,
    nn_param=nn_param,
    bart_param=bart_param
  ) 
  )
  
}


#' Obtain the cross-validated lambda hyparameter for the LASSO
#' @description Obtain the lambda hyparameter for the LASSO using cross-validation 
#' @param folds number of folds
#' @param fmla formula object ex. "E ~ x1+x2"
#' @param tao time point of interest
#' @param data a training data set
#' @return lambda to be used in the glmnet function 
#' @rdname stackBagg-internal



tune_lasso <- function(folds,
                       fmla,
                       tao,
                       data,
                       xnam) {
  
  xnam.factor <- colnames(data[xnam])[sapply(data[xnam], class)=="factor"]
  if(length(xnam.factor)==0){ xnam.factor<- NULL}
  xnam.cont <- xnam[!(xnam %in% xnam.factor)]
  xnam.cont.gam <- xnam.cont[apply(data[xnam.cont],2, function(z) length(unique(z))>3 )]
  
  
  grid.lasso <- glmnetUtils::glmnet(fmla,data=data[!is.na(data$E),],family="binomial")$lambda
  lambda_max <- max(grid.lasso)
  epsilon <- .0001
  K <- 100
  lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = K)), digits = 10)
  lasso_param <- lambdapath
  
  
  pred.test.lasso <- NULL
  id.set <- NULL
  
  test.id_temp=1:nrow(data)
  for (k in 1:folds) {
    
    if (k==folds){
      test.id=test.id_temp
    }else{
      test.id=sample(test.id_temp,floor( nrow(data)/folds ), replace=FALSE) # nrow(data)/folds=100
      test.id_temp <- setdiff(test.id_temp, test.id) # update the test id candidates
    }
    train.id <- setdiff(1:nrow(data), test.id)
    train.set <- data.frame(data[train.id,])
    test.set <- data.frame(data[test.id,])
    n_test.set <- nrow(test.set)
    
    train.set=train.set[!is.na(train.set$E),]
    
    pred.lasso=apply(as.matrix(lasso_param),1,function(x) ML_list$lassofun(data=train.set,testdat=test.set,fmla,xnam,xnam.factor,xnam.cont,x) ) 
    pred.test.lasso=rbind(pred.test.lasso, pred.lasso )
    
    id.set=rbind(id.set,as.matrix(test.set$id))  
  }
  
  pred.test.lasso=pred.test.lasso[order(id.set),]
  data=data[order(data$id),]
  loss.function.lasso=apply(pred.test.lasso,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts,tao))
  
  lasso_param=lasso_param[which.max(loss.function.lasso)]
  return(lasso_param)
}
