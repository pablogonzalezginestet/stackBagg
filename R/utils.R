
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


############### Machine Learning Procedures ############################

#' @description  Predictions based on a library of Machine Learning procedures
#' @param traindata training data set
#' @param testdata validation/test data set
#' @param fmla formula object ex. "E ~ x1+x2"
#' @param tuneparams a list of tune parameters for each machine learning procedure
#' @param i sample selected by bootstrap
#' @return a matrix of predictions where each column is the prediction of each algorithm based on the testdata
#' @rdname EnsBagg-internal


MLprocedures <- function(traindata,testdata,fmla,tuneparams,i){ 
  
  sampledata<- as.data.frame(traindata[i, ])
  
  pred1 <- ML_list$logfun(sampledata,testdata,fmla)
  pred2 <- ML_list$GAMfun(sampledata,testdata,fmla,tuneparams$gam)
  pred3 <- ML_list$lassofun(sampledata,testdata,fmla,tuneparams$lasso)
  pred4 <- ML_list$rffun(sampledata,testdata,fmla,tuneparams$rf)
  pred5 <- ML_list$svmfun(sampledata,testdata, tuneparams$svm)
  pred6 <- ML_list$bartfun(sampledata,testdata, tuneparams$bart)
  pred7 <- ML_list$knnfun(sampledata,testdata, tuneparams$knn)
  pred8 <- ML_list$nn(sampledata,testdata,tuneparams$nn)
  return(cbind(pred1,pred2,pred3,pred4,pred5,pred6,pred7,pred8))
}

#' @description  Library of Machine Learning procedures
#' @importFrom caret dummyVars
#' @return a list of Machine Learning functions
#' @rdname EnsBagg-internal

ML_list <- list(
  
  logfun= function(data,testdata,fmla){
    fit <- stats::glm(fmla, data = data,family = "binomial")
    pred <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
    return(pred)
  },
  
  GAMfun = function(data,testdata,fmla,df) {
    if (length(df)>1){
      
      newterms=c(paste0("s(",xnam.cont.gam, ",df=3)"),attr(terms(fmla[]), "term.labels")[!(attr(terms(fmla[]), "term.labels") %in% xnam.cont.gam)]) 
      newfmla=reformulate(newterms,fmla[[2]])
      fit <- gam(newfmla, data = data, family = 'quasibinomial')
      pred.df3 <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
      
      newterms=c(paste0("s(",xnam.cont.gam, ",df=4)"),attr(terms(fmla[]), "term.labels")[!(attr(terms(fmla[]), "term.labels") %in% xnam.cont.gam)]) 
      newfmla=reformulate(newterms,fmla[[2]])
      fit <- gam(newfmla, data = data, family = 'quasibinomial')
      pred.df4 <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
      return(cbind(pred.df3,pred.df4))
      
    }else{
      
    if (df==3){
      newterms=c(paste0("s(",xnam.cont.gam, ",df=3)"),attr(terms(fmla[]), "term.labels")[!(attr(terms(fmla[]), "term.labels") %in% xnam.cont.gam)]) 
      newfmla=reformulate(newterms,fmla[[2]])
    }else{
      newterms=c(paste0("s(",xnam.cont.gam, ",df=4)"),attr(terms(fmla[]), "term.labels")[!(attr(terms(fmla[]), "term.labels") %in% xnam.cont.gam)]) 
      newfmla=reformulate(newterms,fmla[[2]]) 
    }
    fit <- gam::gam(newfmla, data = data, family = 'quasibinomial')
    pred <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
    return(pred)
    }
    
    }, 
  
  lassofun=function(data,testdata,fmla,lambda){
    fit <- glmnet::glmnet(fmla,data=data,lambda=lambda,family="binomial")
    pred <- predict(fit, newdata = testdata, type = "response", s =lambda)
    return(pred)
  } ,
  
  rffun=function(data,testdata,fmla,grid.rf){
    data=na.omit(cbind(E=as.factor(data$E),data[xnam]))
    fit<- ranger::ranger(fmla,data =data,probability = TRUE, num.trees = grid.rf[1],mtry = grid.rf[2] )
    pred<- predict(fit, data = testdata,type = "response")$predictions[,2]
    return(pred)
  },
  
  svmfun=function(data,testdata,svm.grid){
    Y=na.omit(data$E)
    if(is.null(xnam.factor)){
      X=data[xnam] 
      testdata <- as.data.frame(testdata)[xnam]
    }else{
      X=data[xnam.cont]
      X=cbind(X,predict(dummyVars( ~ ., data =data[xnam.factor], levelsOnly = FALSE), newdata=data[xnam.factor]))
      testdata=as.data.frame(cbind(testdata[xnam.cont],predict(dummyVars( ~ ., data =testdata[xnam.factor], levelsOnly = FALSE), newdata=testdata[xnam.factor])))
    }
    
    X=X[!is.na(data$E),]
    fit.svm <- e1071::svm(y = as.factor(Y), x = X, 
                          type = "C-classification", fitted = FALSE, probability = TRUE, 
                          kernel = "radial", cost =svm.grid[1],gamma=svm.grid[2])
    pred <- attr(predict(fit.svm, newdata = testdata, probability = TRUE), 
                 "prob")[, "1"]
    return(pred)
  } ,
  
  bartfun=function(data,testdata,bart.grid){
    Y=na.omit(data$E)
    X=data[xnam][!is.na(data$E),]
    testdata <- as.data.frame(testdata)[xnam]
    fit <- bartMachine::bartMachine(X,factor(Y, levels = c("1", "0")),num_trees = bart.grid[1],num_burn_in = 250, verbose = FALSE, alpha = 0.95,
                                    beta = 2, k = bart.grid[2], q=bart.grid[3], num_iterations_after_burn_in = 1000)
    pred <- predict(fit,testdata,type="prob" )
    return(pred)
  } ,
  
  knnfun=function(data,testdata,k){
    Y=na.omit(data$E)
    if(is.null(xnam.factor)){
      X=data[xnam] 
      testdata <- as.data.frame(testdata)[xnam]
    }else{
      X=data[xnam.cont]
      X=cbind(X,predict(dummyVars( ~ ., data =data[xnam.factor], levelsOnly = FALSE), newdata=data[xnam.factor]))
      testdata=as.data.frame(cbind(testdata[xnam.cont],predict(dummyVars( ~ ., data =testdata[xnam.factor], levelsOnly = FALSE), newdata=testdata[xnam.factor])))
    }
    X=X[!is.na(data$E),]
    
    fit <- class::knn(train = X, test = testdata, k = k, cl = Y, prob = TRUE)
    pred <- (as.numeric(fit) - 1) * attr(fit, "prob") + (1 - (as.numeric(fit) - 1)) * (1 - attr(fit, "prob"))
    return(pred)
  } ,
  
  
  nnfun=function(traindata,testdata,neurons){
    
    if(is.null(xnam.factor)){
      testdata <- as.data.frame(testdata)[xnam]
      fmla=as.formula(paste("E ~ ", paste(xnam, collapse= "+")))
    }else{
      X=traindata[xnam.cont]
      X=cbind(X,predict(dummyVars( ~ ., data =traindata[xnam.factor], levelsOnly = FALSE), newdata=traindata[xnam.factor]))
      testdata=as.data.frame(cbind(testdata[xnam.cont],predict(dummyVars( ~ ., data =testdata[xnam.factor], levelsOnly = FALSE), newdata=testdata[xnam.factor])))
      colnames(X)=c(paste("x", 1:(dim(X)[2]), sep=""))
      traindata=cbind(E=traindata$E,X)
      colnames(testdata)=c(paste("x", 1:(dim(testdata)[2]), sep=""))
      fmla=as.formula(paste("E ~ ", paste(names(X), collapse= "+")))
    }
    
    traindata=traindata[!is.na(traindata$E),]
    
    fit=neuralnet::neuralnet(fmla,traindata,hidden =as.numeric(neurons) ,threshold = 0.1,act.fct ="logistic" ,linear.output = FALSE,err.fct = "ce", stepmax=1e+08)
    pred=predict(fit, testdata)[,2]
    return(pred)
    
  }
  
)



############### Natively ############################
#' @description  Predictions based on those Machine Learning procedures in the library that allow for weights to be specified as an argument of the R function. No bagging occurs. This group of algorithms is denoted as  Native Weights
#' @param traindata training data set
#' @param testdata validation/test data set
#' @param fmla formula object ex. "E ~ x1+x2"
#' @param tuneparams a list of tune parameters for each machine learning procedure
#' @return a matrix of predictions where each column is the prediction of each algorithm based on the testdata
#' @rdname EnsBagg-internal

MLprocedures_natively <- function(traindata,testdata,fmla,tuneparams){ 
  traindata<- as.data.frame(traindata)
  wts <- traindata$wts.nn 
  
  pred1 <- ML_list_natively$logfun(traindata,testdata,fmla,wts)
  pred2 <- ML_list_natively$GAMfun(traindata,testdata,fmla,tuneparams$gam,wts)
  pred3 <- ML_list_natively$lassofun(traindata,testdata,fmla,tuneparams$lasso,wts)
  pred4 <- ML_list_natively$rffun(traindata,testdata,fmla,tuneparams$rf,wts)
  return(cbind(pred1,pred2,pred3,pred4))
}



ML_list_natively <- list(
  
  
  logfun= function(traindata,testdata,fmla,wts){
    fit <- stats::glm(fmla, data = traindata,family = "binomial", weights=wts)
    pred <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
    return(pred)
  },
  
  
  GAMfun = function(traindata,testdata,fmla,df,wts) {
    if (df==3){
      newterms=c(paste0("s(",xnam.cont.gam, ",df=3)"),attr(terms(fmla[]), "term.labels")[!(attr(terms(fmla[]), "term.labels") %in% xnam.cont.gam)]) 
      newfmla=reformulate(newterms,fmla[[2]])
    }else{
      newterms=c(paste0("s(",xnam.cont.gam, ",df=4)"),attr(terms(fmla[]), "term.labels")[!(attr(terms(fmla[]), "term.labels") %in% xnam.cont.gam)]) 
      newfmla=reformulate(newterms,fmla[[2]]) 
    }
    
    fit <- gam::gam(newfmla, data = traindata, family = 'quasibinomial', weights=wts)
    pred <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
    return(pred)
  } , 
  
  lassofun=function(traindata,testdata,fmla,lambda,wts){
    fit <- glmnet::glmnet(fmla,data=traindata,lambda=lambda,family="binomial", weights=wts)
    pred <- predict(fit, newdata = testdata, type = "response", s =lambda)
    return(pred)
  } ,
  

  rffunw=function(traindata,testdata,fmla,grid.rf,wts){
    data.rf=na.omit(cbind(y=traindata$E,wts=wts,traindata[xnam]))
    fit<- ranger::ranger(y~ .-wts, data =data.rf ,case.weights=data.rf$wts,probability = TRUE,  num.trees =grid.rf[1] ,mtry = grid.rf[2])
    pred<- predict(fit, data = testdata,type = "response")$predictions[,2]
    return(pred)
  }
  
)






