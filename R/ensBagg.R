
#'  Ensemble IPCW Bagging 
#' @description  Main Algorithm
#' @param train.data a data.frame with at least the following variables and column names: 
#' id: unique identifiers of each subject
#' ttilde: event-times (censored), 
#' delta: event indicator at the corresponding value of the variable ttilde. Censored observations must be denoted by the value 0. Main event of interest is denoted by 1. 
#' covariates/features:  these are the covariate that the user potentially want to use in building the preodiction model. 
#' @param test.data a data.frame with the same variables and names that the train.data  
#' @param xnam vector with the names of the covariates to be included in the model
#' @param tao evaluation time point of interest
#' @param weighting CoxPH or CoxBoost 
#' @param folds Number of folds
#' @param tuneparams a list of tune parameters for each machine learning procedure. Name them as gam_param, lasso_param, randomforest_param, svm_param, bart_param, knn_param, nn_param.
#' Default values are the same used for the simulation.
#' @param B number of bootstrap samples
#' @return a list with the predictions of each machine learning algorithm, the average AUC across folds for each of them, the optimal coefficients, an indicator if the optimization procedure has converged and the value of penalization term chosen
#' @importFrom dplyr mutate
#' @import survival
#' @rdname ensbagg
#' @export

ensBagg <- function(train.data,test.data, xnam, tao , weighting , folds , tuneparams=NULL ,B=NULL ){

  
# global parameters
if (missing(B)) {
  B <- 10
}

if (missing(tuneparams)) {
tuneparams <- EnsBagg::parametersSimulation(folds = 5,xnam,train.data)
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

xnam.factor <- colnames(train.data[xnam])[sapply(train.data[xnam], class)=="factor"]
if(length(xnam.factor)==0){ xnam.factor<- NULL}
xnam.cont <- xnam[!(xnam %in% xnam.factor)]
xnam.cont.gam <- xnam.cont[apply(train.data[xnam.cont],2, function(z) length(unique(z))>3 )]


# create binary outcome,E, was created by dichotomizing the time to failure at tao
train.data<- dplyr::mutate(train.data,E=as.factor(ifelse(ttilde < tao & delta==1, 1 , ifelse(ttilde < tao & delta==2 | ttilde>tao, 0, NA))),
                         deltac=ifelse(delta==0,1,0))
test.data<- dplyr::mutate(test.data,E=as.factor(ifelse(ttilde < tao & delta==1, 1 , ifelse(ttilde < tao & delta==2 | ttilde>tao, 0, NA))),
                           deltac=ifelse(delta==0,1,0))

fmla <- as.formula(paste("E ~ ", paste(xnam, collapse= "+")))

if(weighting=="CoxBoost"){
# CoxBoost Weights
#train data
X.boost <- as.matrix(train.data[xnam])
fit.cboost <- peperr::fit.CoxBoost(Surv(train.data$ttilde,train.data$deltac), x=X.boost,cplx=300) 

wts.boost=NULL
for (i in 1:nrow(train.data) ){
  tao_temp <- min(train.data$ttilde[i],tao)
  wts.boost<- c(wts.boost,abs(train.data$deltac[i]-1) / peperr::predictProb(fit.cboost, Surv(train.data$ttilde[i],train.data$deltac[i]), train.data[xnam],tao_temp,complexity = 300)[i])
  }

train.data$wts <- wts.boost
train.data$sum_wts_one <- wts.boost/sum(wts.boost)

#test data
X.boost <- as.matrix(test.data[xnam])
fit.cboost <- peperr::fit.CoxBoost(Surv(test.data$ttilde,test.data$deltac), x=X.boost,cplx=300) 

wts.boost=NULL
for (i in 1:nrow(test.data) ){
  tao_temp <- min(test.data$ttilde[i],tao)
  wts.boost<- c(wts.boost,abs(test.data$deltac[i]-1) / peperr::predictProb(fit.cboost, Surv(test.data$ttilde[i],test.data$deltac[i]), test.data[xnam],tao_temp,complexity = 300)[i])
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
  
  train.data$wts <- wts.coxph
  train.data$sum_wts_one <- wts.coxph/sum(wts.coxph)
  
  #test data set
  cox.C.test<- survival::coxph(fmla.c, data=test.data)
  newdata_zero <- test.data[1,]
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

  test.data$wts <- wts.coxph
  test.data$sum_wts_one <- wts.coxph/sum(wts.coxph)
  
}

train.data <- train.data[c("id","E","wts","sum_wts_one","ttilde","delta",xnam)]
test.data <- test.data[c("id","E","wts","sum_wts_one","ttilde","delta",xnam)]

auc_ipcwBagg <- matrix(NA, nrow = 1 , ncol = A + 1 )
auc_native_weights <- matrix(NA, nrow = 1 , ncol =A_native)

algorithm2<- ipcw_ensbagg(folds=folds, 
                          MLprocedures=MLprocedures,
                          fmla=fmla,
                          tuneparams=tuneparams,
                          B=B,
                          data=train.data ,
                          A=A,
                          xnam=xnam,
                          xnam.factor=xnam.factor,
                          xnam.cont=xnam.cont,
                          xnam.cont.gam=xnam.cont.gam)



prediction_ipcwBagg<- ipcw_genbagg(fmla=fmla,
                                   tuneparams=tuneparams,
                                   MLprocedures=MLprocedures,
                                   traindata = train.data,
                                   testdata = test.data ,
                                   A=A,
                                   xnam=xnam,
                                   xnam.factor=xnam.factor,
                                   xnam.cont=xnam.cont,
                                   xnam.cont.gam=xnam.cont.gam )


prediction_ens_ipcwBagg <- as.matrix(prediction_ipcwBagg) %*% algorithm2$coefficients #combining predictions

auc_ipcwBagg[1,1:A] <- apply(prediction_ipcwBagg,2, function(x) ipcw_auc(T=test.data$ttilde,delta=test.data$delta,marker=crossprod(t(x),1),cause=1,wts=test.data$wts,tao))
auc_ipcwBagg[1,A+1] <- ipcw_auc(T=test.data$ttilde,delta=test.data$delta,marker=prediction_ens_ipcwBagg,cause=1,wts=test.data$wts,tao)

colnames(prediction_ipcwBagg) <- ml_names
colnames(prediction_ens_ipcwBagg) <- c("Ensemble")
colnames(auc_ipcwBagg) <- c(ml_names,"Ensemble")

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
colnames(auc_native_weights) <- c(ml_names[1:A_native])
colnames(prediction_native_weights) <- c(ml_names[1:A_native])

return(list( 
    prediction_ensBagg=cbind(prediction_ipcwBagg,prediction_ens_ipcwBagg),
    auc_ipcwBagg=auc_ipcwBagg,
    prediction_native_weights=prediction_native_weights,
    auc_native_weights=auc_native_weights,
    tuneparams=tuneparams
    ) )

}