#' Tuning parameter selection
#' @description Tuning parameter selection
#' @param gam_param a vector containing degree of freedom 3 and 4
#' @param lasso_param a grid of values for the shrinkage term 
#' @param randomforest_param a two column matrix: first column denotes the num_trees parameter and the second column denotes the mtry parameter.
#' @param knn_param a grid of  positive integers values
#' @param svm_param a three column matrix:  first column denotes the cost parameter, second column the gamma and third column the kernel. kernel=1 denotes "radial" and kernel=2 denotes "linear".
#' @param nn_param a grid of positive integers values for the neurons
#' @param bart_param a three column matrix:  first column denotes the num_tree parameter, second column the k parameter and third column the q parameter.
#' @param folds number of folds
#' @param xnam vector with the names of the covariates to be included in the model
#' @param tao evaluation time point of interest 
#' @param data data set that contains at least id, E , ttilde, delta, wts and covariates
#' @param weighting Procedure to compute the inverse probability of censoring weights. Weighting="CoxPH" and weighting="CoxBoost" model the censoring by the Cox model and CoxBoost model respectively.
#' @return a list with the tune parameters selected using the IPCW AUC loss function.
#' @export



tune_params_ml <- function( gam_param,
                            lasso_param,
                            randomforest_param,
                            knn_param,
                            svm_param,
                            nn_param,
                            bart_param,
                            folds,
                            xnam,
                            tao,
                            data,
                            weighting) {

  xnam.factor <- colnames(data[xnam])[sapply(data[xnam], class)=="factor"]
  if(length(xnam.factor)==0){ xnam.factor<- NULL}
  xnam.cont <- xnam[!(xnam %in% xnam.factor)]
  # continous variables to be applied the smoothing function in the gam must have
  # 10 or more unique values or
  # have to have four values with more than 5%
  xnam.cont.gam <- xnam.cont[apply(data[xnam.cont],2, function(z) length(unique(z))>10 | length(unique(z))<=10 & sum(table(z)/dim(data)[1]>0.05)>=4)]
  
  
  names(data)[1] <- "ttilde"
  names(data)[2] <- "delta"
 

  # create binary outcome,E, was created by dichotomizing the time to failure at tao
  # and we create id: unique identifiers of each subject
  data<- dplyr::mutate(data,id = 1:length(ttilde),E=as.factor(ifelse(ttilde < tao & delta==1, 1 , ifelse(ttilde < tao & delta==2 | ttilde>tao, 0, NA))),
                             deltac=ifelse(delta==0,1,0))
  
  #check that there is no rare categories in each factor variable
  if(!is.null(xnam.factor)) {
    for(s in 1:length(xnam.factor)) {
      level_to_drop <-  table(data[xnam.factor][,s])/dim(data)[1]<.05
      for(q in 1:length(level_to_drop)){
        if(level_to_drop[q]==TRUE){  
          data<- data[!data[xnam.factor][,s]==levels(data[xnam.factor][,s])[q],]
         
        }
      }
      if(any(level_to_drop==TRUE)){
        data[xnam.factor][,s] <- droplevels(data[xnam.factor][,s])
       }
    }
  }
  
  
  
  if(weighting=="CoxBoost"){
    # CoxBoost Weights
    #train data
    if(is.null(xnam.factor)){
      
      X.boost <- as.matrix(data[xnam])
      fit.cboost <- peperr::fit.CoxBoost(Surv(data$ttilde,data$deltac), x=X.boost,cplx=300) 
      wts.boost=NULL
      for (i in 1:nrow(data) ){
        tao_temp <- min(data$ttilde[i],tao)
        wts.boost<- c(wts.boost,as.numeric(!is.na(data$E[i])) / peperr::predictProb(fit.cboost, Surv(data$ttilde[i],data$deltac[i]), data[xnam],tao_temp,complexity = 300)[i])
      }
      
    }else{
      
      X.boost <- data[xnam.cont]
      X.boost <- cbind(X.boost,predict(caret::dummyVars( ~ ., data =data[xnam.factor], fullRank=TRUE), newdata=data[xnam.factor]))
      colnames(X.boost)=c(paste("x", 1:(dim(X.boost)[2]), sep=""))
      data2=as.data.frame(cbind(ttilde=data$ttilde,deltac=data$deltac,X.boost))
      X.boost=as.matrix(data2[,-(1:2)])
      fit.cboost <- peperr::fit.CoxBoost(Surv(data2$ttilde,data2$deltac), x=X.boost,cplx=300) 
      wts.boost=NULL
      xnam2 <-names(data2)[-(1:2)] 
      for (i in 1:nrow(data2) ){
        tao_temp <- min(data2$ttilde[i],tao)
        wts.boost<- c(wts.boost,as.numeric(!is.na(data$E[i])) / peperr::predictProb(fit.cboost, Surv(data2$ttilde[i],data2$deltac[i]), data2[xnam2],tao_temp,complexity = 300)[i])
      }
      
    }
    
    data$wts <- wts.boost
    }
  
  
  if(weighting=="CoxPH"){
    #Cox PH weights
    
    fmla.c <- as.formula(paste("Surv(ttilde,deltac) ~ ", paste(xnam, collapse= "+")))
    #train data set
    cox.C.train<- survival::coxph(fmla.c, data=data)
    newdata_zero <- data[1,]
    
    if(is.null(xnam.factor)){
      
      newdata_zero[xnam] <- 0
      basel_surv <- survival::survfit(cox.C.train, newdata=newdata_zero[xnam])
      beta.cox.train <-  cox.C.train$coef
      
      wts.coxph=NULL
      for (i in 1:nrow(data)){
        newdata <- data[i,]
        newdata_cov <-  data[i,xnam]
        G <- basel_surv$surv ^(exp(sum(newdata_cov * beta.cox.train )))
        wts.coxph <- c(wts.coxph, as.numeric(!is.na(newdata$E)) / (G[ max(which(basel_surv$time <= min(newdata$ttilde, tao) ) ) ]  ) )
      }
      
    }else{
      
      for(s in 1:length(xnam.factor)) {newdata_zero[xnam.factor[s]]=levels(data[xnam.factor][,s])[1]}
      newdata_zero[xnam.cont] <- 0
      basel_surv <- survival::survfit(cox.C.train, newdata=newdata_zero[xnam])
      beta.cox.train <-  cox.C.train$coef
      dummies.factor.train=as.data.frame(predict(caret::dummyVars( ~ ., data =data[xnam.factor],fullRank = TRUE), newdata=data[xnam.factor]))
      beta.cox.train=c(beta.cox.train[xnam.cont],beta.cox.train[!names(beta.cox.train) %in% xnam.cont ])
      
      wts.coxph=NULL
      for (i in 1:nrow(data)){
        newdata <- data[i,]
        newdata_cov = cbind(data[i,xnam.cont],dummies.factor.train[i,])
        G <- basel_surv$surv ^(exp(sum(newdata_cov * beta.cox.train )))
        wts.coxph <- c(wts.coxph, as.numeric(!is.na(newdata$E)) / (G[ max(which(basel_surv$time <= min(newdata$ttilde, tao) ) ) ]  ) )
      }
    }
    
    data$wts <- wts.coxph
    }
  
  
  data <- data[c("id","E","wts","ttilde","delta",xnam)]

  fmla <- as.formula(paste("E ~ ", paste(xnam, collapse= "+")))
  
pred.test.gam=NULL
pred.test.knn=NULL
pred.test.svm=NULL
pred.test.nn=NULL
pred.test.rf=NULL
pred.test.bart=NULL
pred.test.lasso=NULL
id.set=NULL

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
  
  pred.gam=ML_list$GAMfun(data=train.set,
                          testdata=test.set,
                          fmla,
                          xnam,
                          xnam.factor,
                          xnam.cont,
                          xnam.cont.gam,
                          param=gam_param)
  pred.test.gam<-rbind(pred.test.gam, pred.gam )
  
  pred.lasso=apply(as.matrix(lasso_param),1,function(x) ML_list$lassofun(data=train.set,
                                                                         testdata=test.set,
                                                                         fmla,
                                                                         xnam,
                                                                         xnam.factor,
                                                                         xnam.cont,
                                                                         param=x) ) 
  pred.test.lasso=rbind(pred.test.lasso, pred.lasso )
  
  pred.rf=apply(as.matrix(randomforest_param),1,function(x) ML_list$rffun(data=train.set,
                                                                          testdata=test.set,
                                                                          fmla,
                                                                          xnam,
                                                                          xnam.factor,
                                                                          xnam.cont,
                                                                          param=x) ) 
  pred.test.rf=rbind(pred.test.rf, pred.rf ) 
  
  pred.knn=apply(as.matrix(knn_param),1,function(x) ML_list$knnfun(data=train.set,
                                                                   testdata=test.set,
                                                                   xnam,
                                                                   xnam.factor,
                                                                   xnam.cont,
                                                                   param=x)
                                                                   ) 
  pred.test.knn<-rbind(pred.test.knn, pred.knn )
  
  pred.svm=apply(as.matrix(svm_param),1,function(x) ML_list$svmfun(data=train.set,
                                                                   testdata=test.set,
                                                                   xnam,
                                                                   xnam.factor,
                                                                   xnam.cont,
                                                                   param=x)
                                                                   )
  pred.test.svm=rbind(pred.test.svm, pred.svm )
  
  pred.nn=apply(as.matrix(nn_param),1,function(x)  ML_list$nnfun(data=train.set,
                                                              testdata=test.set,
                                                              xnam,
                                                              xnam.factor,
                                                              xnam.cont,
                                                              param=x)
                                                              ) 
  pred.test.nn=rbind(pred.test.nn, pred.nn )
  
  pred.bart=apply(as.matrix(bart_param),1,function(x) ML_list$bartfun(data=train.set,
                                                                      testdata=test.set,
                                                                      xnam,
                                                                      xnam.factor,
                                                                      xnam.cont,
                                                                      param=x)
                                                                      ) 
  pred.test.bart=rbind(pred.test.bart, pred.bart )
  
  
  id.set=rbind(id.set,as.matrix(test.set$id))  
  
  print(k) 
  
}

pred.test.gam=pred.test.gam[order(id.set),]
pred.test.knn=pred.test.knn[order(id.set),]
pred.test.svm=pred.test.svm[order(id.set),]
pred.test.nn=pred.test.nn[order(id.set),]
pred.test.rf=pred.test.rf[order(id.set),]
pred.test.bart=pred.test.bart[order(id.set),]
pred.test.lasso=pred.test.lasso[order(id.set),]
data=data[order(data$id),]

loss.function.gam=apply(pred.test.gam,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts,tao))
loss.function.knn=apply(pred.test.knn,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts,tao))
loss.function.svm=apply(pred.test.svm,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts,tao))
loss.function.nn=apply(pred.test.nn,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts,tao))
loss.function.rf=apply(pred.test.rf,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts,tao))
loss.function.bart=apply(pred.test.bart,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts,tao))
loss.function.lasso=apply(pred.test.lasso,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts,tao))

tuneparams=list(
  gam_param=gam_param[which.max(loss.function.gam)],
  knn_param=knn_param[which.max(loss.function.knn)],
  svm_param=svm_param[which.max(loss.function.svm),],
  nn_param=nn_param[which.max(loss.function.nn)],
  randomforest_param=randomforest_param[which.max(loss.function.rf),],
  bart_param=bart_param[which.max(loss.function.bart),],
  lasso_param=lasso_param[which.max(loss.function.lasso)]
)

  names(tuneparams$lasso_param)=c("lambda")
  names(tuneparams$knn_param)=c("k")
  names(tuneparams$gam_param)=c("df")
  names(tuneparams$nn_param)=c("neurons")
  names(tuneparams$randomforest_param)=c("num_trees","mtry")
  names(tuneparams$svm_param)=c("cost","gamma","kernel")
  names(tuneparams$bart_param)=c("num_tree","k","q")

  return(tuneparams)
  
}











