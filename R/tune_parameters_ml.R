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
                            data ) {
  
data<- dplyr::mutate(data,id = 1:length(ttilde),E=as.factor(ifelse(ttilde < tao & delta==1, 1 , ifelse(ttilde < tao & delta==2 | ttilde>tao, 0, NA))) )

xnam.factor <- colnames(data[xnam])[sapply(data[xnam], class)=="factor"]
if(length(xnam.factor)==0){ xnam.factor<- NULL}
xnam.cont <- xnam[!(xnam %in% xnam.factor)]
xnam.cont.gam <- xnam.cont[apply(data[xnam.cont],2, function(z) length(unique(z))>10 | length(unique(z))<=10 & sum(table(z)/dim(data)[1]>0.05)>=4)]
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
  gam=gam_param[which.max(loss.function.gam)],
  knn=knn_param[which.max(loss.function.knn)],
  svm=svm_param[which.max(loss.function.svm),],
  nn=nn_param[which.max(loss.function.nn)],
  rf=randomforest_param[which.max(loss.function.rf),],
  bart=bart_param[which.max(loss.function.bart),],
  lasso=lasso_param[which.max(loss.function.lasso)]
)

  names(tuneparams$lasso)=c("lambda")
  names(tuneparams$knn)=c("k")
  names(tuneparams$gam)=c("df")
  names(tuneparams$nn)=c("neurons")
  names(tuneparams$rf)=c("num_trees","mtry")
  names(tuneparams$svm)=c("cost","gamma")
  names(tuneparams$bart)=c("num_tree","k","q")

  return(tuneparams)
  
}











