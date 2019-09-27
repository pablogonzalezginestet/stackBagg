
#'  Ensemble IPCW Bagging 
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
#' @param tuneparams a list of tune parameters for each machine learning procedure. Name them as gam_param, lasso_param, randomforest_param, svm_param, bart_param, knn_param, nn_param.
#' Default values are the same used for the simulation.
#' @param B number of bootstrap samples
#' @return a list with the predictions of each machine learning algorithm, the average AUC across folds for each of them, the optimal coefficients, an indicator if the optimization procedure has converged and the value of penalization term chosen
#' @import survival
#' @rdname ensbagg
#' @export

ensBagg <- function(train.data,test.data, xnam, tao , weighting , folds , tuneparams=NULL ,B=NULL,discard=NULL ){

  
if (missing(B)) {
  B <- 10
}

xnam.factor <- colnames(train.data[xnam])[sapply(train.data[xnam], class)=="factor"]
if(length(xnam.factor)==0){ xnam.factor<- NULL}
xnam.cont <- xnam[!(xnam %in% xnam.factor)]
xnam.cont.gam <- xnam.cont[apply(train.data[xnam.cont],2, function(z) length(unique(z))>3 )]

# checking if the data was provided in the right form
#if(names(train.data)[1]!="id" | names(train.data)[1]!="ID" | names(train.data)[1]!="Id" ){
#  stop("column id is missing. Column id must be the first column. Check appropiate format in the help file")
#}
#if(names(test.data)[1]!="id" | names(test.data)[1]!="ID" | names(test.data)[1]!="Id" ){
 # stop("column id is missing. Column id must be the first column. Check appropiate format in the help file")
#}
# rename second and third column
#names(train.data)[1] <- "id"
names(train.data)[1] <- "ttilde"
names(train.data)[2] <- "delta"
#names(test.data)[1] <- "id"
names(test.data)[1] <- "ttilde"
names(test.data)[2] <- "delta"


# create binary outcome,E, was created by dichotomizing the time to failure at tao
# and we create id: unique identifiers of each subject
train.data<- dplyr::mutate(train.data,id = 1:length(ttilde),E=as.factor(ifelse(ttilde < tao & delta==1, 1 , ifelse(ttilde < tao & delta==2 | ttilde>tao, 0, NA))),
                         deltac=ifelse(delta==0,1,0))
test.data<- dplyr::mutate(test.data,id = 1:length(ttilde),E=as.factor(ifelse(ttilde < tao & delta==1, 1 , ifelse(ttilde < tao & delta==2 | ttilde>tao, 0, NA))),
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



auc_ipcwBagg <- matrix(NA, nrow = 1 , ncol = A + 1 )
auc_native_weights <- matrix(NA, nrow = 1 , ncol =A_native)

algorithm2<- ipcw_ensbagg( folds=folds, 
                          MLprocedures=MLprocedures,
                          fmla=fmla,
                          tuneparams=tuneparams,
                          tao,
                          B=B,
                          data=train.data ,
                          A=A,
                          xnam=xnam,
                          xnam.factor=xnam.factor,
                          xnam.cont=xnam.cont,
                          xnam.cont.gam=xnam.cont.gam )


 
prediction_ipcwBagg<- ipcw_genbagg( fmla=fmla,
                                   tuneparams=tuneparams,
                                   MLprocedures=MLprocedures,
                                   traindata = train.data,
                                   testdata = test.data ,
                                   A=A,
                                   B=B,
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

# Survival Methods

#survival
# cause specific Cox regression
fmla.cox <- as.formula(paste("Hist(ttilde,delta) ~ ", paste(xnam, collapse= "+")))
coxfit <- riskRegression::CSC(formula = fmla.cox,data = train.data) 
predcoxfit <- riskRegression::predictRisk(coxfit, newdata = test.data[xnam], cause = 1, times = tao)
auc_survival1 <- ipcw_auc(T=test.data$ttilde,delta=test.data$delta,marker=predcoxfit,cause=1,wts=test.data$wts,tao)

#cox boost 
X.coxboost=as.matrix(train.data[xnam])
newX.coxboost <- as.matrix(as.data.frame(test.data)[xnam])
fitCoxboost <- CoxBoost::CoxBoost(time=train.data$ttilde,status=train.data$delta, x=X.coxboost,stepno=300,penalty=100)
predCoxboost <- predict(fitCoxboost,newdata=newX.coxboost,times=tao,type="CIF")
auc_survival2 <- ipcw_auc(T=test.data$ttilde,delta=test.data$delta,marker=predCoxboost,cause=1,wts=test.data$wts,tao)

#random forest
fmla.rf <- as.formula(paste("Surv(ttilde,delta) ~ ", paste(xnam, collapse= "+")))
fitSurvRf <- randomForestSRC::rfsrc(fmla.rf ,data=train.data ,nsplit = 3,ntree = 100)
predSurvRf <- randomForestSRC::predict.rfsrc(fitSurvRf, test.data[xnam])
time.index.cif.rf <- data.table::last(which(predSurvRf$time.interest<25))
cifRf <- predSurvRf$cif[,time.index.cif.rf,1]
auc_survival3 <- ipcw_auc(T=test.data$ttilde,delta=test.data$delta,marker=cifRf,cause=1,wts=test.data$wts,tao)

auc_survival<- cbind(auc_survival1,auc_survival2,auc_survival3)
colnames(auc_survival) <- c("CoxPH","CoxBoost","Random Forest")
prediction_survival <- cbind(predcoxfit,predCoxboost,cifRf)
colnames(prediction_survival) <- c("CoxPH","CoxBoost","Random Forest")

if(!missing(discard)){
  dat.discard.train <- train.data[!is.na(train.data$E),]
  dat.discard.test <- test.data[!is.na(test.data$E),]
  
  pred1 <- ML_list$logfun(dat.discard.train,dat.discard.test,fmla, xnam,xnam.factor,xnam.cont)
  pred2 <- ML_list$GAMfun(dat.discard.train,dat.discard.test,fmla,xnam,xnam.factor,xnam.cont,xnam.cont.gam,tuneparams$gam_param)
  pred3 <- ML_list$lassofun(dat.discard.train,dat.discard.test,fmla,xnam,xnam.factor,xnam.cont,tuneparams$lasso_param)
  pred4 <- ML_list$rffun(dat.discard.train,dat.discard.test,fmla,xnam,xnam.factor,xnam.cont,tuneparams$randomforest_param)
  pred5 <- ML_list$svmfun(dat.discard.train,dat.discard.test,xnam,xnam.factor,xnam.cont,tuneparams$svm_param)
  pred6 <- ML_list$bartfun(dat.discard.train,dat.discard.test,xnam,xnam.factor,xnam.cont, tuneparams$bart_param)
  pred7 <- ML_list$knnfun(dat.discard.train,dat.discard.test,xnam,xnam.factor,xnam.cont, tuneparams$knn_param)
  pred8 <- ML_list$nn(dat.discard.train,dat.discard.test,xnam,xnam.factor,xnam.cont,tuneparams$nn_param)
  prediction_discard<- as.matrix(cbind(pred1,pred2,pred3,pred4,pred5,pred6,pred7,pred8))
  colnames(prediction_discard) <- ml_names
  auc_discard <- apply(prediction_discard,2, function(x) cvAUC::AUC(predictions = x,labels =dat.discard.test$E) )
  
  return(list( 
    prediction_ensBagg=cbind(prediction_ipcwBagg,prediction_ens_ipcwBagg),
    prediction_native_weights=prediction_native_weights,
    prediction_survival=prediction_survival,
    prediction_discard=prediction_discard,
    optimal_coefficients=algorithm2$coefficients,
    algorithm2$convergence_indicator,
    algorithm2$penalization_term,
    auc_ipcwBagg=auc_ipcwBagg,
    auc_native_weights=auc_native_weights,
    auc_survival=auc_survival,
    auc_discard=auc_discard,
    tuneparams=tuneparams
  ) )
}

if(missing(discard)){

return(list( 
    prediction_ensBagg=cbind(prediction_ipcwBagg,prediction_ens_ipcwBagg),
    prediction_native_weights=prediction_native_weights,
    prediction_survival=prediction_survival,
    optimal_coefficients=algorithm2$coefficients,
    algorithm2$convergence_indicator,
    algorithm2$penalization_term,
    auc_ipcwBagg=auc_ipcwBagg,
    auc_native_weights=auc_native_weights,
    auc_survival=auc_survival,
    tuneparams=tuneparams
    ) )
}

}