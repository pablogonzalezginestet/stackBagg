
#' Internal stackBagg functions
#' 
#' @description Internal stackBagg helper functions
#'
#' @details These functions are not intended for use by users.
#'
#' @name stackBagg-internal

NULL

#' @rdname stackBagg-internal
#' @param coef_init  starting values for the coefficients
#' @param lambda penalization term. It is a positive scalar.
#' @param data A data frame that contains at least: ttilde, delta, wts
#' @param Z a matrix that contains the predictions. Each column represents a single marker.
#' @param tao time point of interest 
#' @return a vector with the optimal AUC value and the optimal coefficient  


optimun_auc_coef = function(coef_init,lambda,data,Z,tao){
  optimal_coef <- optim(par = coef_init, fn = risk_auc,lambda=lambda ,Z=Z, data=data,tao=tao,method = "BFGS",control=list(maxit=10000))
  AUC_coef_opt<- ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(Z),optimal_coef$par),cause=1,wts=data$wts,tao)
  return(c(AUC_coef_opt,optimal_coef$convergence,optimal_coef$par))
}


#' @description  Compute the risk of missclassifying an individual using as a marker a single prediction or weighted linear combination of several predictions (1-AUC) 
#' @param par  a vector of coefficients/weights. Its length must be equal to the number of predictions included in Z 
#' @param lambda penalization term. It is a positive scalar.
#' @param Z a matrix that contains the predictions. Each column represents a single marker.
#' @param data A data frame  that constains at least: ttilde= time to event, delta=event type, wts= IPC weights
#' @param tao time point of interest 
#' @return 1-AUC
#' @rdname stackBagg-internal


risk_auc<- function(par,lambda,Z,data,tao){
  par <- par/sum(par) 
  marker_pred <- crossprod(t(Z),par)
  Risk_AUC=1-ipcw_auc(T=data$ttilde,delta=data$delta,marker=marker_pred,cause=1,wts=data$wts,tao)+lambda*sum(abs(par))
  return(Risk_AUC)
}


############### Machine Learning Procedures ############################

#' Library of Machine learning procedures
#' @description  Library of Machine Learning procedures
#' @return a character vector with all prediction algorithms  supported by stackBagg 
#' @export

  algorithms <- function(){
  return(c("ens.glm","ens.gam","ens.lasso","ens.randomForest","ens.svm","ens.bartMachine","ens.knn","ens.nn"))
}



#' @description  Predictions based on a library of Machine Learning procedures
#' @param traindata training data set
#' @param testdata validation/test data set
#' @param fmla formula object ex. "E ~ x1+x2"
#' @param tuneparams a list of tune parameters for each machine learning procedure
#' @param i sample selected by bootstrap
#' @return a matrix of predictions where each column is the prediction of each algorithm based on the testdata
#' @rdname stackBagg-internal

#MLprocedures <- function(traindata,testdata,fmla, tuneparams,i){ 
  
  
 # sampledata<- as.data.frame(traindata[i, ])
  
  #pred1 <- ML_list$logfun(sampledata,testdata,fmla)
  #pred2 <- ML_list$GAMfun(sampledata,testdata,fmla, tuneparams$gam_param)
  #pred3 <- ML_list$lassofun(sampledata,testdata,fmla, tuneparams$lasso_param)
  #pred4 <- ML_list$rffun(sampledata,testdata,fmla, tuneparams$randomforest_param)
  #pred5 <- ML_list$svmfun(sampledata,testdata, tuneparams$svm_param)
  #pred6 <- ML_list$bartfun(sampledata,testdata, tuneparams$bart_param)
  #pred7 <- ML_list$knnfun(sampledata,testdata, tuneparams$knn_param)
  #pred8 <- ML_list$nn(sampledata,testdata, tuneparams$nn_param)
 # return(cbind(pred1,pred2,pred3,pred4,pred5,pred6,pred7,pred8))
#}



MLprocedures <- function(traindata,
                         testdata,
                         fmla,
                         xnam,
                         xnam.factor,
                         xnam.cont,
                         xnam.cont.gam,
                         tuneparams,
                         ens.library,
                         i){ 
  
  
  sampledata<- as.data.frame(traindata[i, ])
  pred <- NULL
  if("ens.glm" %in% ens.library){
    pred_temp <- ML_list$logfun(sampledata,testdata,fmla, xnam,xnam.factor,xnam.cont)
    pred <- cbind(pred_temp,pred)
  }
  if("ens.gam" %in% ens.library){
    pred_temp <- ML_list$GAMfun(sampledata,testdata,fmla,xnam,xnam.factor,xnam.cont,xnam.cont.gam,tuneparams$gam_param)
    pred <- cbind(pred_temp,pred)
  }
  if("ens.lasso" %in% ens.library){
    pred_temp <- ML_list$lassofun(sampledata,testdata,fmla,xnam,xnam.factor,xnam.cont,tuneparams$lasso_param)
    pred <- cbind(pred_temp,pred)
  }
  if("ens.randomForest" %in% ens.library){
    pred_temp <- ML_list$rffun(sampledata,testdata,fmla,xnam,xnam.factor,xnam.cont,tuneparams$randomforest_param)
    pred <- cbind(pred_temp,pred)
  }
  if("ens.svm" %in% ens.library){
    pred_temp <- ML_list$svmfun(sampledata,testdata,xnam,xnam.factor,xnam.cont,tuneparams$svm_param)
    pred <- cbind(pred_temp,pred)
  }
  
  if("ens.bartMachine" %in% ens.library){
    pred_temp <-  ML_list$bartfun(sampledata,testdata,xnam,xnam.factor,xnam.cont, tuneparams$bart_param)
    pred <- cbind(pred_temp,pred)
  }
  
  
  if("ens.knn" %in% ens.library){
    pred_temp <-  ML_list$knnfun(sampledata,testdata,xnam,xnam.factor,xnam.cont, tuneparams$knn_param)
    pred <- cbind(pred_temp,pred)
  }
  
  if("ens.nn" %in% ens.library){
    pred_temp <- ML_list$nn(sampledata,testdata,xnam,xnam.factor,xnam.cont,tuneparams$nn_param)
    pred <- cbind(pred_temp,pred)
  }
  
  return(pred)  
  
}

  
  
#' @description  Library of Machine Learning procedures
#' @return a list of Machine Learning functions
#' @rdname stackBagg-internal

ML_list <- list(
  
  logfun= function(data,testdata,fmla, xnam,xnam.factor,xnam.cont){
    fit <- stats::glm(fmla, data = data,family = "binomial")
    pred <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
    return(pred)
  },
 
  GAMfun = function(data,
                    testdata,
                    fmla,
                    xnam,
                    xnam.factor,
                    xnam.cont,
                    xnam.cont.gam,
                    param) {
  
        if (length(param)>1){
      
      newterms=c(paste0("s(",xnam.cont.gam, ",df=3)"),xnam[!xnam %in% xnam.cont.gam]) 
      newfmla=stats::reformulate(newterms,fmla[[2]])
      fit1 <- gam::gam(newfmla, data = data, family = 'quasibinomial')
      pred.df3 <- predict(fit1, newdata=testdata, type = "response", na.action=na.omit)
     
      newterms=c(paste0("s(",xnam.cont.gam, ",df=4)"),xnam[!xnam %in% xnam.cont.gam]) 
      newfmla=stats::reformulate(newterms,fmla[[2]])
      fit2 <- gam::gam(newfmla, data = data, family = 'quasibinomial')
      pred.df4 <- predict(fit2, newdata=testdata, type = "response", na.action=na.omit)
      pred <- cbind(pred.df3,pred.df4)
      return(pred)
      
    }else if (param==3){
      newterms=c(paste0("s(",xnam.cont.gam, ",df=3)"),xnam[!xnam %in% xnam.cont.gam]) 
      newfmla=stats::reformulate(newterms,fmla[[2]])
      fit <- gam::gam(newfmla, data = data, family = 'quasibinomial')
      pred <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
      return(pred)
      }else {
        
      newterms=c(paste0("s(",xnam.cont.gam, ",df=4)"),xnam[!xnam %in% xnam.cont.gam])
      newfmla=stats::reformulate(newterms,fmla[[2]]) 
      fit <- gam::gam(newfmla, data = data, family = 'quasibinomial')
      pred <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
      return(pred)
      }
   
    }, 
  
  lassofun=function(data,testdata,fmla,xnam,xnam.factor,xnam.cont,param){

    fit <- glmnetUtils::glmnet(fmla,data=data,lambda=param,family="binomial")
    pred <- predict(fit, newdata = testdata, type = "response", s =param)
    return(pred)
  } ,
  
  rffun=function(data,testdata,fmla,xnam,xnam.factor,xnam.cont,param){
  
    data=na.omit(cbind(E=as.factor(data$E),data[xnam]))
    fit<- ranger::ranger(fmla,data =data,probability = TRUE, num.trees = param[1],mtry = param[2] )
    pred<- predict(fit, data = testdata,type = "response")$predictions[,2]
    return(pred)
  },
  
  svmfun=function(data,testdata,xnam,xnam.factor,xnam.cont,param){

    Y=na.omit(data$E)
    if(is.null(xnam.factor)){
      X=data[xnam] 
      testdata <- as.data.frame(testdata)[xnam]
    }else{
      X=data[xnam.cont]
      X=cbind(X,predict(caret::dummyVars( ~ ., data =data[xnam.factor], fullRank=TRUE), newdata=data[xnam.factor]))
      testdata=as.data.frame(cbind(testdata[xnam.cont],predict(caret::dummyVars( ~ ., data =testdata[xnam.factor],  fullRank=TRUE), newdata=testdata[xnam.factor])))
    }
    
    X=X[!is.na(data$E),]
    if(param[3]==1){
    fit.svm <- e1071::svm(y = as.factor(Y), x = X, 
                          type = "C-classification", fitted = FALSE, probability = TRUE, 
                          kernel = "radial" , cost =param[1],gamma=param[2])
    }else{
      fit.svm <- e1071::svm(y = as.factor(Y), x = X, 
                            type = "C-classification", fitted = FALSE, probability = TRUE, 
                            kernel = "linear" , cost =param[1])
    }
    pred <- attr(predict(fit.svm, newdata = testdata, probability = TRUE), 
                 "prob")[, "1"]
    return(pred)
  } ,
  
  bartfun=function(data,testdata,xnam,xnam.factor,xnam.cont,param){
 
   Y=na.omit(data$E)
    X=data[xnam][!is.na(data$E),]
    testdata <- as.data.frame(testdata)[xnam]
    fit <- bartMachine::bartMachine(X,factor(Y, levels = c("1", "0")),num_trees = param[1],num_burn_in = 250, verbose = FALSE, alpha = 0.95,
                                    beta = 2, k = param[2], q=param[3], num_iterations_after_burn_in = 1000)
    pred <- predict(fit,testdata,type="prob" )
    return(pred)
  } ,
  
  knnfun=function(data,testdata,xnam,xnam.factor,xnam.cont,param){

   Y=na.omit(data$E)
    if(is.null(xnam.factor)){
      X=data[xnam] 
      testdata <- as.data.frame(testdata)[xnam]
    }else{
      X=data[xnam.cont]
      X=cbind(X,predict(caret::dummyVars( ~ ., data =data[xnam.factor], fullRank=TRUE), newdata=data[xnam.factor]))
      testdata=as.data.frame(cbind(testdata[xnam.cont],predict(caret::dummyVars( ~ ., data =testdata[xnam.factor],  fullRank=TRUE), newdata=testdata[xnam.factor])))
    }
    X=X[!is.na(data$E),]
    
    fit <- class::knn(train = X, test = testdata, k = param, cl = Y, prob = TRUE)
    pred <- (as.numeric(fit) - 1) * attr(fit, "prob") + (1 - (as.numeric(fit) - 1)) * (1 - attr(fit, "prob"))
    return(pred)
  } ,
  
  
  nnfun=function(data,testdata,xnam,xnam.factor,xnam.cont,param){
  
   
    if(is.null(xnam.factor)){
      testdata <- as.data.frame(testdata)[xnam]
      fmla=as.formula(paste("E ~ ", paste(xnam, collapse= "+")))
    }else{
      X=data[xnam.cont]
      X=cbind(X,predict(caret::dummyVars( ~ ., data =data[xnam.factor], fullRank=TRUE), newdata=data[xnam.factor]))
      testdata=as.data.frame(cbind(testdata[xnam.cont],predict(caret::dummyVars( ~ ., data =testdata[xnam.factor],  fullRank=TRUE), newdata=testdata[xnam.factor])))
      colnames(X)=c(paste("x", 1:(dim(X)[2]), sep=""))
      data=cbind(E=data$E,X)
      colnames(testdata)=c(paste("x", 1:(dim(testdata)[2]), sep=""))
      fmla=as.formula(paste("E ~ ", paste(names(X), collapse= "+")))
    }
    
    traindata=data[!is.na(data$E),]
    
    fit=neuralnet::neuralnet(fmla,traindata,hidden =as.numeric(param) ,threshold = 0.1,act.fct ="logistic" ,linear.output = FALSE,err.fct = "ce", stepmax=1e+08)
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
#' @rdname stackBagg-internal



MLprocedures_natively <- function(
         traindata,
         testdata,
         fmla,
         xnam,
         xnam.factor,
         xnam.cont,
         xnam.cont.gam,
         tuneparams
         ){ 
  
  wts <- traindata$wts
  pred <- NULL
  if("ens.glm" %in% ens.library){
  pred_temp <- ML_list_natively$logfun(traindata,testdata,fmla, xnam,xnam.factor,xnam.cont,wts)
  pred <- cbind(pred_temp,pred)
  }

  if("ens.gam" %in% ens.library){
  pred_temp <- ML_list_natively$GAMfun(traindata,testdata,fmla,xnam,xnam.factor,xnam.cont,xnam.cont.gam,tuneparams$gam_param,wts)
  pred <- cbind(pred_temp,pred)
  }
  
  if("ens.lasso" %in% ens.library){
  pred_temp <- ML_list_natively$lassofun(traindata,testdata,fmla,xnam,xnam.factor,xnam.cont,tuneparams$lasso_param,wts)
  pred <- cbind(pred_temp,pred)
  }
  
  if("ens.randomForest" %in% ens.library){
  pred_temp <- ML_list_natively$rffun(traindata,testdata,fmla,xnam,xnam.factor,xnam.cont,tuneparams$randomforest_param)
  pred <- cbind(pred_temp,pred)
  }
  
  return(pred)
}


#' @description  Library of Machine Learning procedures that allows for weights
#' @return a list of Machine Learning functions
#' @rdname stackBagg-internal
#' 
ML_list_natively <- list(
  
  logfun= function(traindata,testdata,fmla, xnam,xnam.factor,xnam.cont,wts){
    fit <- stats::glm(fmla, data = traindata,family = "binomial", weights=wts)
    pred <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
    return(pred)
  },
  
  
  GAMfun = function(traindata,
                    testdata,
                    fmla,
                    xnam,
                    xnam.factor,
                    xnam.cont,
                    xnam.cont.gam,
                    param,
                    wts) {
    if (length(param)>1){
      
      newterms=c(paste0("gam::s(",xnam.cont.gam, ",df=3)"),xnam[!xnam %in% xnam.cont.gam]) 
      newfmla=stats::reformulate(newterms,fmla[[2]])
      
      fit <- gam::gam(newfmla, data = traindata, family = 'quasibinomial',weights=wts)
      pred.df3 <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
      
      newterms=c(paste0("gam::s(",xnam.cont.gam, ",df=4)"),xnam[!xnam %in% xnam.cont.gam]) 
      newfmla=stats::reformulate(newterms,fmla[[2]])
      fit <- gam::gam(newfmla, data = traindata, family = 'quasibinomial',weights=wts)
      pred.df4 <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
      
      return(cbind(pred.df3,pred.df4))
      
    }else{
      
      if (param==3){
        newterms=c(paste0("gam::s(",xnam.cont.gam, ",df=3)"),xnam[!xnam %in% xnam.cont.gam]) 
        newfmla=stats::reformulate(newterms,fmla[[2]])
      }else{
        newterms=c(paste0("gam::s(",xnam.cont.gam, ",df=4)"),xnam[!xnam %in% xnam.cont.gam])
        newfmla=stats::reformulate(newterms,fmla[[2]]) 
      }
      fit <- gam::gam(newfmla, data = traindata, family = 'quasibinomial', weights=wts)
      pred <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)
      
      return(pred)
    }
    
  }, 
  
  lassofun=function(traindata,testdata,fmla,xnam,xnam.factor,xnam.cont,param,wts){
    fit <- glmnetUtils::glmnet(fmla,data=traindata,lambda=param,family="binomial", weights=wts)
    pred <- predict(fit, newdata = testdata, type = "response", s =param)
    return(pred)
  } ,
  
  rffun=function(traindata,testdata,fmla,xnam,xnam.factor,xnam.cont,param){
    data=na.omit(cbind(E=as.factor(traindata$E), wts=traindata$wts,traindata[xnam]))
    fit<- ranger::ranger(fmla,data =data,probability = TRUE, case.weights=data$wts , num.trees = param[1],mtry = param[2] )
    pred<- predict(fit, data = testdata,type = "response")$predictions[,2]
    return(pred)
  }
  
)


#########    Tuning Parameter    #################
#' Grid of values for the Real Data Application: InfCareHIV Register
#' @description  A grid of values for hyperparameters used in the Real Data Application: InfCareHIV Register. This grid of values isan argument in the tuning parameter function tune_parameter_ml.R 
#' @param xnam a vector with the covariates names considered in the modeling
#' @param data a training data set
#' @param tao time point of interest
#' @return a list with a grid of values for each hyperparameter
#'  gam_param a vector containing degree of freedom 3 and 4
#'  lasso_param a grid of values for the shrinkage term lambda
#'  randomforest_param a two column matrix: first column denotes the num_trees parameter and the second column denotes the mtry parameter.
#'  knn_param a grid of  positive integers values
#'  svm_param a three column matrix:  first column denotes the cost parameter, second column the gamma and third column the kernel. kernel=1 denotes "radial" and kernel=2 denotes "linear".
#' nn_param a grid of positive integers values for the neurons
#'  bart_param a three column matrix:  first column denotes the num_tree parameter, second column the k parameter and third column the q parameter.
#' @export
#' @rdname stackBagg-internal


grid_parametersDataHIV <- function(xnam,data,tao){
  
  data<- dplyr::mutate(data,E=as.factor(ifelse(ttilde < tao & delta==1, 1 , ifelse(ttilde < tao & delta==2 | ttilde>tao, 0, NA))) )
  
  fmla <- as.formula(paste("E ~ ", paste(xnam, collapse= "+")))
  
  gam_param <- c(3,4)
  
  grid.lasso <-glmnetUtils::glmnet(fmla,data=data[!is.na(data$E),],family="binomial")$lambda
  lambda_max <- max(grid.lasso)
  epsilon <- .0001
  K <- 100
  lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = K)), digits = 10)
  lasso_param <- lambdapath
  
  
  num.trees <-  c(50,100,250,500,750,1000)
  mtry <- c(1,floor(sqrt(ncol(data[xnam]))),floor(ncol(data[xnam])/3),floor(ncol(data[xnam])/3),seq(floor(ncol(data[xnam])/3)+1,19,by=2)) 
  randomforest_param <- cbind(rep(num.trees,times=rep(length(mtry),length(num.trees))),rep(mtry,length(num.trees)))
  
  knn_param <- seq(1,50,by=1)
  
  cost <- c(10^3, 10^2, 10, 1)
  gamma <- c(10^(-5), 10^(-4), 10^(-3), 10^(-2), 10^(-1))
  svm_param <- rbind(cbind(cost=rep(cost, times=rep(length(gamma),length(cost))),gamma=rep(gamma,length(cost)),kernel=1 ), cbind(cost=c(10,1,0.1,.01),gamma=NA,kernel=2) )
  
  nn_param <- c(1,2,3,4,5)
  
  num_tree <-  c(50, 200)
  k <-  c(1,2,3,5)
  q <- c(0.9,0.99,0.75)
  grid.temp <- cbind( rep(num_tree,times=rep(length(k),length(num_tree))), rep(k,length(num_tree))  )
  bart_param <-  cbind(rep(grid.temp[,1],times=rep(length(q),length(grid.temp[,1]))),rep(grid.temp[,2],times=rep(length(q),length(grid.temp[,2]))), 
                    rep(q,length(grid.temp[,1]))  )
  
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
