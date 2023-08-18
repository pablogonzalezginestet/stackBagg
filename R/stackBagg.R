
#' Stacked IPCW Bagging 
#' @description  Main Algorithm
#' @param train.data a data.frame with at least the following variables: 
#' event-times (censored) in the first column,
#' event indicator in the second column and 
#' covariates/features that the user potentially want to use in building the preodiction model.
#' Censored observations must be denoted by the value 0. Main event of interest is denoted by 1.
#' @param test.data a data.frame with the same variables and names that the train.data  
#' @param xnam vector with the names of the covariates to be included in the model
#' @param tao evaluation time point of interest
#' @param weighting Procedure to compute the inverse probability of censoring weights. Weighting="CoxPH" and weighting="CoxBoost" model the censoring by the Cox model and CoxBoost model respectively.
#' @param folds Number of folds
#' @param ens.library character vector indicating the prediction algorithms to be consider in the analyisis. The prediction algorithms supported by this package are: "ens.glm","ens.gam","ens.lasso","ens.randomForest","ens.svm","ens.bartMachine","ens.knn","ens.nn"). See the function ensBagg::ens.all.algorithms(). 
#' @param tuneparams a list of tune parameters for each machine learning procedure. Name them as gam_param, lasso_param, randomforest_param, svm_param, bart_param, knn_param, nn_param.
#' Default values are the same used for the simulation.
#' @param B number of bootstrap samples
#' @return a list with the predictions of each machine learning algorithm, the average AUC across folds for each of them, the optimal coefficients,the weights ,an indicator if the optimization procedure has converged and the value of penalization term chosen
#' @import survival
#' @import gam
#' @import Matrix
#' @rdname stackBagg
#' @export

stackBagg <- function(train.data,test.data, xnam, tao , weighting , folds ,ens.library, tuneparams=NULL ,B=NULL ){

 all.library <- stackBagg::algorithms()

if (missing(B)) {
  B <- 10
}

xnam.factor <- colnames(train.data[xnam])[sapply(train.data[xnam], class)=="factor"]
if(length(xnam.factor)==0){ xnam.factor<- NULL}
xnam.cont <- xnam[!(xnam %in% xnam.factor)]
# continous variables to be applied the smoothing function in the gam must have
# 10 or more unique values or
# have to have four values with more than 5%

xnam.cont.gam <- xnam.cont[apply(train.data[xnam.cont],2, function(z) length(unique(z))>10 | length(unique(z))<=10 & sum(table(z)/dim(train.data)[1]>0.05)>=4)]


# checking if the data was provided in the right form

names(train.data)[1] <- "ttilde"
names(train.data)[2] <- "delta"

names(test.data)[1] <- "ttilde"
names(test.data)[2] <- "delta"


# create binary outcome,E, was created by dichotomizing the time to failure at tao
# and we create id: unique identifiers of each subject
train.data<- dplyr::mutate(train.data,id = 1:length(ttilde),E=as.factor(ifelse(ttilde < tao & delta==1, 1 , ifelse(ttilde < tao & delta==2 | ttilde>tao, 0, NA))),
                         deltac=ifelse(delta==0,1,0))
test.data<- dplyr::mutate(test.data,id = 1:length(ttilde),E=as.factor(ifelse(ttilde < tao & delta==1, 1 , ifelse(ttilde < tao & delta==2 | ttilde>tao, 0, NA))),
                           deltac=ifelse(delta==0,1,0))

#check that there is no rare categories in each factor variable
if(!is.null(xnam.factor)) {
  for(s in 1:length(xnam.factor)) {
    level_to_drop <-  table(train.data[xnam.factor][,s])/dim(train.data)[1]<.05
    for(q in 1:length(level_to_drop)){
      if(level_to_drop[q]==TRUE){  
        train.data<- train.data[!train.data[xnam.factor][,s]==levels(train.data[xnam.factor][,s])[q],]
        test.data<- test.data[!test.data[xnam.factor][,s]==levels(test.data[xnam.factor][,s])[q],]
      }
    }
    if(any(level_to_drop==TRUE)){
      train.data[xnam.factor][,s] <- droplevels(train.data[xnam.factor][,s])
      test.data[xnam.factor][,s]=droplevels(test.data[xnam.factor][,s])}
  }
}



fmla <- as.formula(paste("E ~ ", paste(xnam, collapse= "+")))

if(weighting=="CoxBoost"){
# CoxBoost Weights
#train data
if(is.null(xnam.factor)){
  
  X.boost <- as.matrix(train.data[xnam])
  fit.cboost <- peperr::fit.CoxBoost(Surv(train.data$ttilde,train.data$deltac), x=X.boost,cplx=300) 
  wts.boost=NULL
  for (i in 1:nrow(train.data) ){
    tao_temp <- min(train.data$ttilde[i],tao)
    wts.boost<- c(wts.boost,as.numeric(!is.na(train.data$E[i])) / peperr::predictProb(fit.cboost, Surv(train.data$ttilde[i],train.data$deltac[i]), train.data[xnam],tao_temp,complexity = 300)[i])
  }
  
  }else{
    
  X.boost <- train.data[xnam.cont]
  X.boost <- cbind(X.boost,predict(caret::dummyVars( ~ ., data =train.data[xnam.factor], fullRank=TRUE), newdata=train.data[xnam.factor]))
  colnames(X.boost)=c(paste("x", 1:(dim(X.boost)[2]), sep=""))
  train.data2=as.data.frame(cbind(ttilde=train.data$ttilde,deltac=train.data$deltac,X.boost))
  X.boost=as.matrix(train.data2[,-(1:2)])
  fit.cboost <- peperr::fit.CoxBoost(Surv(train.data2$ttilde,train.data2$deltac), x=X.boost,cplx=300) 
  wts.boost=NULL
  xnam2 <-names(train.data2)[-(1:2)] 
  for (i in 1:nrow(train.data2) ){
    tao_temp <- min(train.data2$ttilde[i],tao)
    wts.boost<- c(wts.boost,as.numeric(!is.na(train.data$E[i])) / peperr::predictProb(fit.cboost, Surv(train.data2$ttilde[i],train.data2$deltac[i]), train.data2[xnam2],tao_temp,complexity = 300)[i])
  }
  
  }

train.data$wts <- wts.boost
train.data$sum_wts_one <- wts.boost/sum(wts.boost)

#test data
if(is.null(xnam.factor)){
X.boost <- as.matrix(test.data[xnam])
fit.cboost <- peperr::fit.CoxBoost(Surv(test.data$ttilde,test.data$deltac), x=X.boost,cplx=300) 

wts.boost=NULL
for (i in 1:nrow(test.data) ){
  tao_temp <- min(test.data$ttilde[i],tao)
  wts.boost<- c(wts.boost,as.numeric(!is.na(test.data$E[i])) / peperr::predictProb(fit.cboost, Surv(test.data$ttilde[i],test.data$deltac[i]), test.data[xnam],tao_temp,complexity = 300)[i])
  }

} else {
  X.boost <- test.data[xnam.cont]
  X.boost <- cbind(X.boost,predict(caret::dummyVars( ~ ., data =test.data[xnam.factor], fullRank=TRUE), newdata=test.data[xnam.factor]))
  colnames(X.boost)=c(paste("x", 1:(dim(X.boost)[2]), sep=""))
  test.data2=as.data.frame(cbind(ttilde=test.data$ttilde,deltac=test.data$deltac,X.boost))
  X.boost=as.matrix(test.data2[,-(1:2)])
  fit.cboost <- peperr::fit.CoxBoost(Surv(test.data2$ttilde,test.data2$deltac), x=X.boost,cplx=300) 
  wts.boost=NULL
  xnam2 <-names(test.data2)[-(1:2)] 
  for (i in 1:nrow(test.data2) ){
    tao_temp <- min(test.data2$ttilde[i],tao)
    wts.boost<- c(wts.boost,as.numeric(!is.na(test.data$E[i])) / peperr::predictProb(fit.cboost, Surv(test.data2$ttilde[i],test.data2$deltac[i]), test.data2[xnam2],tao_temp,complexity = 300)[i])
  }
  
}


test.data$wts <- wts.boost
test.data$sum_wts_one <- wts.boost/sum(wts.boost)

}

if(weighting=="CoxPH"){
  #Cox PH weights
  
  fmla.c <- as.formula(paste("Surv(ttilde,deltac) ~ ", paste(xnam, collapse= "+")))
  #train data set
  cox.C.train<- survival::coxph(fmla.c, data=train.data)
  newdata_zero <- train.data[1,]
  
  if(is.null(xnam.factor)){
    
    newdata_zero[xnam] <- 0
    basel_surv <- survival::survfit(cox.C.train, newdata=newdata_zero[xnam])
    beta.cox.train <-  cox.C.train$coef
    
    wts.coxph=NULL
    for (i in 1:nrow(train.data)){
      newdata <- train.data[i,]
      newdata_cov <-  train.data[i,xnam]
      G <- basel_surv$surv ^(exp(sum(newdata_cov * beta.cox.train )))
      wts.coxph <- c(wts.coxph, as.numeric(!is.na(newdata$E)) / (G[ max(which(basel_surv$time <= min(newdata$ttilde, tao) ) ) ]  ) )
    }
    
  }else{
    
  for(s in 1:length(xnam.factor)) {newdata_zero[xnam.factor[s]]=levels(train.data[xnam.factor][,s])[1]}
  newdata_zero[xnam.cont] <- 0
  basel_surv <- survival::survfit(cox.C.train, newdata=newdata_zero[xnam])
  beta.cox.train <-  cox.C.train$coef
  dummies.factor.train=as.data.frame(predict(caret::dummyVars( ~ ., data =train.data[xnam.factor],fullRank = TRUE), newdata=train.data[xnam.factor]))
  beta.cox.train=c(beta.cox.train[xnam.cont],beta.cox.train[!names(beta.cox.train) %in% xnam.cont ])
  
  wts.coxph=NULL
  for (i in 1:nrow(train.data)){
    newdata <- train.data[i,]
    newdata_cov = cbind(train.data[i,xnam.cont],dummies.factor.train[i,])
    G <- basel_surv$surv ^(exp(sum(newdata_cov * beta.cox.train )))
    wts.coxph <- c(wts.coxph, as.numeric(!is.na(newdata$E)) / (G[ max(which(basel_surv$time <= min(newdata$ttilde, tao) ) ) ]  ) )
  }
  }
  
  train.data$wts <- wts.coxph
  train.data$sum_wts_one <- wts.coxph/sum(wts.coxph)
  
  #test data set
  cox.C.test<- survival::coxph(fmla.c, data=test.data)
  newdata_zero <- test.data[1,]
  
  if(is.null(xnam.factor)){
    
    newdata_zero[xnam]=0
    basel_surv <- survival::survfit(cox.C.test, newdata=newdata_zero[xnam])
    beta.cox.test <-  cox.C.test$coef
    
    wts.coxph=NULL
    for (i in 1:nrow(test.data)){
      newdata <- test.data[i,]
      newdata_cov <-  test.data[i,xnam]
      G <- basel_surv$surv ^(exp(sum(newdata_cov * beta.cox.test )))
      wts.coxph <- c(wts.coxph, as.numeric(!is.na(newdata$E)) / (G[ max(which(basel_surv$time <= min(newdata$ttilde, tao) ) ) ]  ) )
    }
  }else{
    
  for(s in 1:length(xnam.factor)) {newdata_zero[xnam.factor[s]]=levels(test.data[xnam.factor][,s])[1]}
  newdata_zero[xnam.cont] <- 0
  basel_surv <- survival::survfit(cox.C.test, newdata=newdata_zero[xnam])
  beta.cox.test <-  cox.C.test$coef
  dummies.factor.test=as.data.frame(predict(caret::dummyVars( ~ ., data =test.data[xnam.factor],fullRank = TRUE), newdata=test.data[xnam.factor]))
  beta.cox.test=c(beta.cox.test[xnam.cont],beta.cox.test[!names(beta.cox.test) %in% xnam.cont ])
  
  wts.coxph=NULL
  for (i in 1:nrow(test.data)){
    newdata <- test.data[i,]
    newdata_cov = cbind(test.data[i,xnam.cont],dummies.factor.test[i,])
    G <- basel_surv$surv ^(exp(sum(newdata_cov * beta.cox.test )))
    wts.coxph <- c(wts.coxph, as.numeric(!is.na(newdata$E)) / (G[ max(which(basel_surv$time <= min(newdata$ttilde, tao) ) ) ]  ) )
  }
}
  test.data$wts <- wts.coxph
  test.data$sum_wts_one <- wts.coxph/sum(wts.coxph)
  
}


train.data <- train.data[c("id","E","wts","sum_wts_one","ttilde","delta",xnam)]
test.data <- test.data[c("id","E","wts","sum_wts_one","ttilde","delta",xnam)]


if (missing(tuneparams)) {
  tuneparams <- stackBagg::parametersSimulation(folds = 5,xnam,train.data,tao)
}  

all.natively<- c("ens.glm","ens.gam","ens.lasso","ens.randomForest")

if("ens.gam" %in% ens.library & length(tuneparams$gam)>1){
  A<- length(ens.library)+1 # if we consider gam with df=3 and df=4
  A_native <- sum(all.natively %in% ens.library)+1
  #ml_names<- c("LogisticReg","GAM.3","GAM.4","LASSO","Random Forest","SVM","BART","k-NN","Neural Network")
  ml_names<- all.library[all.library %in% ens.library]
  ml_names <- c(ml_names[1],"ens.gam.3","ens.gam.4",ml_names[-(1:2)])
  ml_names_natively <- all.natively[all.natively %in% ens.library]
  ml_names_natively <- c(ml_names_natively[1], "ens.gam.3","ens.gam.4",ml_names_natively[-(1:2)])
  }else{
  A<- length(ens.library)
  A_native <- sum(all.natively %in% ens.library)
  #ml_names<- c("LogisticReg","GAM","LASSO","Random Forest","SVM","BART","k-NN","Neural Network")
  ml_names<- all.library[all.library %in% ens.library]
  ml_names_natively <- all.natively[all.natively %in% ens.library]
  }



auc_ipcwBagg <- matrix(NA, nrow = 1 , ncol = A + 1 )
auc_ipcwBagg_training <- matrix(NA, nrow = 1 , ncol = A + 1 )
auc_native_weights <- matrix(NA, nrow = 1 , ncol =A_native)

algorithm2<- ipcw_ensbagg( folds=folds, 
                          MLprocedures=MLprocedures,
                          fmla=fmla,
                          tuneparams=tuneparams,
                          tao,
                          B=B,
                          A,
                          data=train.data ,
                          xnam=xnam,
                          xnam.factor=xnam.factor,
                          xnam.cont=xnam.cont,
                          xnam.cont.gam=xnam.cont.gam,
                          ens.library=ens.library)


 
prediction_ipcwBagg<- ipcw_genbagg( fmla=fmla,
                                   tuneparams=tuneparams,
                                   MLprocedures=MLprocedures,
                                   traindata = train.data,
                                   testdata = test.data ,
                                   B=B,
                                   A,
                                   xnam=xnam,
                                   xnam.factor=xnam.factor,
                                   xnam.cont=xnam.cont,
                                   xnam.cont.gam=xnam.cont.gam,
                                   ens.library)


prediction_ens_ipcwBagg <- as.matrix(prediction_ipcwBagg) %*% algorithm2$coefficients #combining predictions

auc_ipcwBagg[1,1:A] <- apply(prediction_ipcwBagg,2, function(x) ipcw_auc(T=test.data$ttilde,delta=test.data$delta,marker=crossprod(t(x),1),cause=1,wts=test.data$wts,tao))
auc_ipcwBagg[1,A+1] <- ipcw_auc(T=test.data$ttilde,delta=test.data$delta,marker=prediction_ens_ipcwBagg,cause=1,wts=test.data$wts,tao)

colnames(prediction_ipcwBagg) <- ml_names
colnames(prediction_ens_ipcwBagg) <- c("Stack")
colnames(auc_ipcwBagg) <- c(ml_names,"Stack")

# auc in the training set
prediction_ipcwBagg_training<- ipcw_genbagg( fmla=fmla,
                                             tuneparams=tuneparams,
                                             MLprocedures=MLprocedures,
                                             traindata = train.data,
                                             testdata = train.data ,
                                             B=B,
                                             A,
                                             xnam=xnam,
                                             xnam.factor=xnam.factor,
                                             xnam.cont=xnam.cont,
                                             xnam.cont.gam=xnam.cont.gam,
                                             ens.library)

prediction_ens_ipcwBagg_training <- as.matrix(prediction_ipcwBagg_training) %*% algorithm2$coefficients #combining predictions

auc_ipcwBagg_training[1,1:A] <- apply(prediction_ipcwBagg_training,2, function(x) ipcw_auc(T=train.data$ttilde,delta=train.data$delta,marker=crossprod(t(x),1),cause=1,wts=train.data$wts,tao))
auc_ipcwBagg_training[1,A+1] <- ipcw_auc(T=train.data$ttilde,delta=train.data$delta,marker=prediction_ens_ipcwBagg_training,cause=1,wts=train.data$wts,tao)

colnames(auc_ipcwBagg_training) <- c(ml_names,"Stack")


# Native Weights
prediction_native_weights <- MLprocedures_natively(   traindata = train.data,
                                                      testdata = test.data ,
                                                      fmla=fmla,
                                                      xnam,
                                                      xnam.factor,
                                                      xnam.cont,
                                                      xnam.cont.gam,
                                                      tuneparams )
                                     
  
auc_native_weights[1,1:A_native] <- apply(prediction_native_weights,2, function(x) ipcw_auc(T=test.data$ttilde,delta=test.data$delta,marker=crossprod(t(x),1),cause=1,wts=test.data$wts,tao))
auc_native_weights <- as.matrix(auc_native_weights)
colnames(auc_native_weights) <- ml_names_natively
colnames(prediction_native_weights) <- ml_names_natively

# Survival Methods

#survival
# cause specific Cox regression
fmla.cox <- as.formula(paste("Hist(ttilde,delta) ~ ", paste(xnam, collapse= "+")))
coxfit <- riskRegression::CSC(formula = fmla.cox,data = train.data,cause=1,surv.type="surv") 
predcoxfit <- tryCatch( riskRegression::predictRisk(coxfit, newdata = test.data[xnam], cause = 1, times = tao), error=function(e) { e ; return(NA) } )
if (all(!is.nan(predcoxfit))) { 
auc_survival1 <- stackBagg::ipcw_auc(T=test.data$ttilde,delta=test.data$delta,marker=predcoxfit,cause=1,wts=test.data$wts,tao)
}else{
  auc_survival1 <- NA
}
#cox boost 

#if(is.null(xnam.factor)){
#  X.coxboost=as.matrix(train.data[xnam])
#  newX.coxboost <- as.matrix(as.data.frame(test.data)[xnam])
#  fitCoxboost <- CoxBoost::CoxBoost(time=train.data$ttilde,status=train.data$delta, x=X.coxboost,stepno=300,penalty=100)
#  pred.coxboost <- predict(fitCoxboost,newdata=newX.coxboost,times=tao,type="CIF")

#}else{
#X.cboost.train <- train.data[xnam.cont]
#X.cboost.train<- cbind(X.cboost.train,predict(caret::dummyVars( ~ ., data =train.data[xnam.factor], fullRank = TRUE), newdata=train.data[xnam.factor]))
#colnames(X.cboost.train)=c(paste("x", 1:(dim(X.cboost.train)[2]), sep=""))
#train.data2=cbind(ttilde=train.data$ttilde,delta=train.data$delta,X.cboost.train)

#X.cboost.test <- test.data[xnam.cont]
#X.cboost.test <- cbind(X.cboost.test,predict(caret::dummyVars( ~ ., data =test.data[xnam.factor],fullRank = TRUE), newdata=test.data[xnam.factor]))
#colnames(X.cboost.test)=c(paste("x", 1:(dim(X.cboost.test)[2]), sep=""))
#test.data2=cbind(ttilde=test.data$ttilde,delta=test.data$delta,X.cboost.test)

#X.coxboost=as.matrix(train.data2[c(-1,-2)])
#newX.coxboost <- as.matrix(as.data.frame(test.data2)[c(-1,-2)])
#fit.coxboost <-  CoxBoost::CoxBoost(time=train.data2$ttilde,status=train.data2$delta, x=X.coxboost,stepno=300,penalty=100)
#pred.coxboost <- predict(fit.coxboost,newdata=newX.coxboost,times=tao,type="CIF")
#}

#auc_survival2 <- ipcw_auc(T=test.data$ttilde,delta=test.data$delta,marker=pred.coxboost,cause=1,wts=test.data$wts,tao)

#random forest
fmla.rf <- as.formula(paste("Surv(ttilde,delta) ~ ", paste(xnam, collapse= "+")))
fitSurvRf <- randomForestSRC::rfsrc(fmla.rf ,data=train.data ,nsplit = 3,ntree = 100)
predSurvRf <- randomForestSRC::predict.rfsrc(fitSurvRf, test.data[xnam])
time.index.cif.rf <- data.table::last(which(predSurvRf$time.interest<25))
cifRf <- predSurvRf$cif[,time.index.cif.rf,1]
auc_survival3 <- ipcw_auc(T=test.data$ttilde,delta=test.data$delta,marker=cifRf,cause=1,wts=test.data$wts,tao)

#auc_survival<- cbind(auc_survival1,auc_survival2,auc_survival3)
#colnames(auc_survival) <- c("CoxPH","CoxBoost","Random Forest")
#prediction_survival <- cbind(predcoxfit,pred.coxboost,cifRf)
#colnames(prediction_survival) <- c("CoxPH","CoxBoost","Random Forest")

auc_survival<- cbind(auc_survival1,auc_survival3)
colnames(auc_survival) <- c("CoxPH","Random Forest")
prediction_survival <- cbind(predcoxfit,cifRf)
colnames(prediction_survival) <- c("CoxPH","Random Forest")

return(list( 
  library=ml_names,
    prediction_ensBagg=cbind(prediction_ipcwBagg,prediction_ens_ipcwBagg),
    prediction_native_weights=prediction_native_weights,
    prediction_survival=prediction_survival,
    optimal_coefficients=algorithm2$coefficients,
    convergence=algorithm2$convergence_indicator,
    penalization_term=algorithm2$penalization_term,
    auc_ipcwBagg_training=auc_ipcwBagg_training,
    auc_ipcwBagg=auc_ipcwBagg,
    auc_native_weights=auc_native_weights,
    auc_survival=auc_survival,
    wts_train=train.data$wts,
    wts_test=test.data$wts,
    tuneparams=tuneparams
    ) )

}