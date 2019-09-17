
############################# Algorithm 2 Ensemble IPCW Bagging ##############################################
#' @description  Obtain predictions
#' @param folds Number of folds
#' @param MLprocedures \link{MLprocedures}
#' @param fmla formula object ex. "E ~ x1+x2"
#' @param tuneparams a list of tune parameters for each machine learning procedure
#' @param B number of bootstrap samples
#' @param data a training data set
#' @param A a number of machine learning algorithms in the library
#' @return a list with the predictions of each machine learning algorithm and the average AUC across folds for each of them



ipcw_ensbagg <- function(folds, MLprocedures, fmla, tuneparams , B=NULL, data,A) {
  result <- vector("list", A)
  AUC.train <- vector("list", A)
  result_id <- vector("list")
  test.id_temp=1:nrow(data)
  for (k in 1:folds) {
    if (k==folds){
      test.id=test.id_temp
    }else{
      test.id=sample(test.id_temp,floor( nrow(data)/folds ), replace=FALSE)
      test.id_temp <- setdiff(test.id_temp, test.id) # update the test id candidates
    }
    train.id <- setdiff(1:nrow(data), test.id)
    train.set <- data.frame(data[train.id,])
    test.set <- data.frame(data[test.id,])
    n_test.set <- nrow(test.set)
    
    #boot
    b <-boot::boot(data=train.set, statistic=MLprocedures, R=B, fmla=fmla,tuneparams=tuneparams,
                   testdata=test.set, weights = train.set$wts)
    
    d<- apply(b$t,1,function(x) split(x, rep(seq(A), each = n_test.set)))
    D.all <- list()
    d.all <- NULL
    for(s in 1:A){
      for(j in 1:B){
        d.all <- rbind(d.all,unlist(d[[j]][s],use.names = FALSE)) #given a ml procedure i go over the bootstrap for that given ml procedure
      }
      D.all <- colMeans(d.all,na.rm=T)
      result[[s]]<-c(result[[s]], D.all)
      AUC.train[[s]] <- c(AUC.train[[s]], AUC_function(T=test.set$ttilde,delta=test.set$delta,marker=D.all,cause=1,wts=test.set$wts.nn,tao))
      d.all <- NULL
    }
    result_id[[k]]<- test.set$id
    
    # update progress bar
    print(k)
  }
  
  matrix = do.call(cbind, result)
  prediction=as.data.frame(cbind(cbind(unlist(result_id)),matrix))
  prediction=prediction[order(prediction[,1]),]
  names(prediction) = c("id",paste("Set", 1:(ncol(prediction)-1), sep=""))
  AUC.train <- sapply(AUC.train, function(x) sum(x)/folds)
  data <- data[order(data$id),] # sort the data by id since the predictions are sorted by id too 
  
  coef_init <- rep(1/ncol(prediction),ncol(prediction)) #initial values
  penal_grid=c(.01,.1,.5,1,5,10,15,25,50,100) # grid of values of the penalization term considered in the optimization problem
  auc_coef <- matrix(NA,length(penal_grid),ncol(prediction)+2) # a matrix that store the AUC value at the optimum coefficients, if it has converged, the penalization term selected and the optimum coefficients
  for(i in 1:length(penal_grid)){
    auc_coef[i,] <- optimun_auc_coef(penal_grid[i],data,prediction)
  }
  
  coef_opt <- auc_coef[which.max(auc_coef[,1]),][-(1:2)]
  coef_opt_normalized <- coef_opt/sum(coef_opt) # optimal coefficients normalized 
  convergence_indicator <- auc_coef[which.max(auc_coef[,1]),][2]
  penal_chosen <- penal_grid[which.max(auc_coef[,1])]
  
  
  return(list(prediction,AUC.train,coef_opt_normalized,convergence_indicator,penal_chosen))
  
}



############################# Algorithm 1 General IPCW Bagging ##############################################
#' @description  Obtain predictions
#' @param fmla formula object ex. "E ~ x1+x2"
#' @param tuneparams a list of tune parameters for each machine learning procedure
#' @param MLprocedures \link{MLprocedures}
#' @param traindata a training data set
#' @param testdata a test data set 
#' @return a list with the predictions of each machine learning algorithm and the average AUC across folds for each of them
#' @rdname wBrierScore
#' @export

pred_function <- function(fmla,tuneparams,MLprocedures,traindata,testdata) {
  
  result <- vector("list", A)
  n_testdata <- nrow(testdata)
  b <-boot(data=traindata, statistic=MLprocedures, R=B, fmla=formula,tune.params=tune.params,
           testdata=testdata, weights = traindata$wts)
  
  d<- apply(b$t,1,function(x) split(x, rep(seq(A), each = n_testdata)))
  D.all <- list()
  d.all <- NULL
  for(s in 1:A){
    for(j in 1:B){
      d.all <- rbind(d.all,unlist(d[[j]][s],use.names = FALSE)) #given a ml procedure i go over the bootstrap for that given ml procedure
    }
    D.all <- colMeans(d.all,na.rm=T)
    result[[s]]<-c(result[[s]], D.all)
    d.all <- NULL
  }
  
  prediction = do.call(cbind, result)
  return(prediction)
}