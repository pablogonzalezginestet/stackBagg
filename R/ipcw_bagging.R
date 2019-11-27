
#' Algorithm 2:  Procedure to obtain optimally the coefficients to be used in Algorithm 1
#' @description  Obtain predictions
#' @param folds Number of folds
#' @param MLprocedures \link{MLprocedures}
#' @param fmla formula object ex. "E ~ x1+x2"
#' @param tuneparams a list of tune parameters for each machine learning procedure
#' @param tao time point of interest
#' @param B number of bootstrap samples
#' @param data a training data set
#' @param xnam all covariates in the model
#' @param xnam.factor categorical variables include in the model
#' @param xnam.cont continous variables include in the model
#' @param xnam.cont.gam continous variables to be included in the smoothing operator gam::s(,df)
#' @param ens.library  algorithms in the library
#' @return a list with the predictions of each machine learning algorithm (id, predictions), the average AUC across folds for each of them, the optimal coefficients, an indicator if the optimization procedure has converged and the value of penalization term chosen
#' @rdname stackBagg-internal


ipcw_ensbagg <- function(folds,
                         MLprocedures,
                         fmla,
                         tuneparams ,
                         tao,
                         B=NULL,
                         A,
                         data ,
                         xnam,
                         xnam.factor,
                         xnam.cont,
                         xnam.cont.gam,
                         ens.library) {
  
  
  result <- vector("list", A)
  AUC.train <- vector("list", A)
  result_id <- vector("list")
  test.id_temp=1:nrow(data)
  # create progress bar
  pb= txtProgressBar(min = 0, max = folds, style = 3, char=":)")
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
   b <-boot::boot(data=train.set,
                  testdata=test.set,
                  fmla=fmla,
                  xnam=xnam,
                  xnam.factor=xnam.factor,
                  xnam.cont=xnam.cont,
                  xnam.cont.gam=xnam.cont.gam,
                  tuneparams=tuneparams,
                  ens.library=ens.library,
                  statistic=MLprocedures,
                  R=B,
                  weights = train.set$sum_wts_one )
    
    d<- apply(b$t,1,function(x) split(x, rep(seq(A), each = n_test.set)))
    D.all <- list()
    d.all <- NULL
    for(s in 1:A){
      for(j in 1:B){
        d.all <- rbind(d.all,unlist(d[[j]][s],use.names = FALSE)) #given a ml procedure i go over the bootstrap for that given ml procedure
      }
      D.all <- colMeans(d.all,na.rm=T)
      result[[s]]<-c(result[[s]], D.all)
      AUC.train[[s]] <- c(AUC.train[[s]], ipcw_auc(T=test.set$ttilde,delta=test.set$delta,marker=D.all,cause=1,wts=test.set$wts,tao))
      d.all <- NULL
    }
    result_id[[k]]<- test.set$id
    
    # update progress bar
    setTxtProgressBar(pb, k)
  }
  
  matrix = do.call(cbind, result)
  prediction=as.data.frame(cbind(cbind(unlist(result_id)),matrix))
  prediction=prediction[order(prediction[,1]),]
  names(prediction) = c("id",paste("Set", 1:(ncol(prediction)-1), sep=""))
  AUC.train <- sapply(AUC.train, function(x) sum(x)/folds)
  data <- data[order(data$id),] # sort the data by id since the predictions are sorted by id too 
  
  coef_init <- rep(1/A,A) #initial values
  penal_grid=c(.01,.1,.5,1,5,10,15,25,50,100) # grid of values of the penalization term considered in the optimization problem
  auc_coef <- matrix(NA,length(penal_grid),ncol(prediction[,-1])+2) # a matrix that store the AUC value at the optimum coefficients, if it has converged, the penalization term selected and the optimum coefficients
  for(i in 1:length(penal_grid)){
    auc_coef[i,] <- optimun_auc_coef(coef_init,penal_grid[i],data,prediction[,-1],tao)
  }
  
  coef_opt <- auc_coef[which.max(auc_coef[,1]),][-(1:2)]
  coef_opt_normalized <- coef_opt/sum(coef_opt) # optimal coefficients normalized 
  convergence_indicator <- auc_coef[which.max(auc_coef[,1]),][2]
  penal_chosen <- penal_grid[which.max(auc_coef[,1])]
  
  
  
  return(list(prediction=prediction,auc=AUC.train,coefficients=coef_opt_normalized,convergence_indicator=convergence_indicator,penalization_term=penal_chosen))
  
}



#' Algorithm 1: Stacked IPCW Bagging 

#' @description  Obtain predictions
#' @importFrom boot boot
#' @param fmla formula object ex. "E ~ x1+x2"
#' @param tuneparams a list of tune parameters for each machine learning procedure
#' @param MLprocedures \link{MLprocedures}
#' @param traindata a training data set
#' @param testdata a test data set 
#' @param B number of bootstrap samples
#' @param xnam all covariates in the model
#' @param xnam.factor categorical variables include in the model
#' @param xnam.cont continous variables include in the model
#' @param xnam.cont.gam continous variables to be included in the smoothing operator gam::s(,df=)
#' @param ens.library  algorithms in the library
#' @return a matrix with the predictions on the test data set of each machine learning algorithm considered in \link{MLprocedures}
#' @rdname stackBagg-internal


ipcw_genbagg <- function(fmla,
                         tuneparams,
                         MLprocedures,
                         traindata,
                         testdata,
                         B,
                         A,
                         xnam,
                         xnam.factor,
                         xnam.cont,
                         xnam.cont.gam,
                         ens.library) {
  
  result <- vector("list", A)
  n_testdata <- nrow(testdata)
  
  b <-boot::boot(data=traindata,
                 testdata=testdata,
                 fmla=fmla,
                 xnam=xnam,
                 xnam.factor=xnam.factor,
                 xnam.cont=xnam.cont,
                 xnam.cont.gam=xnam.cont.gam,
                 tuneparams=tuneparams,
                 ens.library=ens.library,
                 statistic=MLprocedures,
                 R=B,
                 weights = traindata$sum_wts_one )
  

  
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
