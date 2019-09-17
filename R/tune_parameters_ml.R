

#setwd("/Users/pgonzalezginestet/Documents/KarolinskaInstitutet/1st_Project/Code/Hebbe/My_version_withErinsuggestion/Scenario1")
rm(list = ls())
setwd("Z:/Project 1/Data")
load("hiv.train.RData")
source("functions_packages.R")
#load("hiv.test.RData")

data=hiv.train
B=10 #boostraps
folds=5 
Tstar=tao=730

xnam=c("age","sex","immigrant","ethnicity","infection_route","count_languages","condition.at.supression","count_drug.at.supression",names(hiv.train)[15:31] ) #i took away birht country since it has too many factors
xnam.factor=colnames(data)[sapply(data, class)=="factor"]
xnam.cont=xnam[!(xnam %in% xnam.factor)]
xnam.cont.gam=xnam.cont[sapply(apply(hiv.train[xnam.cont],2,unique),function(z) length(z)>3)]
fmla <- as.formula(paste("E ~ ", paste(xnam, collapse= "+")))

grid.gam=c(3,4)
grid.knn=seq(1,50,by=1)
#Radial basis function
cost=c(10^3, 10^2, 10, 1)
gamma=c(10^(-5), 10^(-4), 10^(-3), 10^(-2), 10^(-1))
grid.svm=rbind(cbind(cost=rep(cost, times=rep(length(gamma),length(cost))),gamma=rep(gamma,length(cost)),kernel=1 ), cbind(cost,gamma=NA,kernel=2) )

grid.nn=c(1,2,3,4,5) #,floor( length(xnam)/2 ))   #seq(1,20,by=1)
num.trees = c(50,100,250,500,750,1000)
mtry=c(1,floor(sqrt(ncol(data[xnam]))),floor(ncol(data[xnam])/3),floor(ncol(data[xnam])/3),seq(floor(ncol(data[xnam])/3)+1,19,by=2)) 
grid.rf=cbind(rep(num.trees,times=rep(length(mtry),length(num.trees))),rep(mtry,length(num.trees)))

#default used in the function bartMachineCV
num_tree = c(50, 200)
k = c(1,2,3,5)
#nu = c(3,10) #not used for classification 
q=c(0.9,0.99,0.75)
grid.temp=cbind( rep(num_tree,times=rep(length(k),length(num_tree))), rep(k,length(num_tree))  )

grid.bart= cbind(rep(grid.temp[,1],times=rep(length(q),length(grid.temp[,1]))),rep(grid.temp[,2],times=rep(length(q),length(grid.temp[,2]))), 
                 rep(q,length(grid.temp[,1]))  )
grid.lasso=glmnet(fmla,data=data[!is.na(data$E),],family="binomial")$lambda
lambda_max <- max(grid.lasso)
epsilon <- .0001
K <- 100
lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), length.out = K)), digits = 10)
grid.lasso=lambdapath


tune_params_ml <- function( gam_param,
                            lasso_param,
                            randomforest_param,
                            knn_param,
                            svm_param,
                            nn_param,
                            bart_param,
                            folds,
                            fmla,
                            tao,
                            data ) {

pred.test.gam=NULL
pred.test.knn=NULL
pred.test.svm=NULL
pred.test.nn=NULL
pred.test.rf=NULL
pred.test.bart=NULL
pred.test.lasso=NULL
id.set=NULL

test.id_temp=1:nrow(data)
#data$E=as.numeric(data$E)-1

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
  
  pred.gam=ML_list$GAMfun(data=train.set,testdata=test.set,fmla,gam_param)
  pred.test.gam<-rbind(pred.test.gam, pred.gam )
  
  pred.lasso=apply(as.matrix(lasso_param),1,function(x) ML_list$lassofun(data=train.set,testdata=test.set,fmla,x) ) 
  pred.test.lasso=rbind(pred.test.lasso, pred.lasso )
  
  pred.rf=apply(as.matrix(randomforest_param),1,function(x) ML_list$rffun(data=train.set,testdata=test.set,fmla,x) ) 
  pred.test.rf=rbind(pred.test.rf, pred.rf ) 
  
  pred.knn=apply(as.matrix(knn_param),1,function(x) ML_list$knnfun(data=train.set,testdata=test.set, x) ) 
  pred.test.knn<-rbind(pred.test.knn, pred.knn )
  
  pred.svm=apply(as.matrix(svm_param),1,function(x) ML_list$svmfun(data=train.set,testdata=test.set, x) )
  pred.test.svm=rbind(pred.test.svm, pred.svm )
  
  pred.nn=apply(as.matrix(nn_param),1,function(x)  ML_list$nn(data=train.set,testdata=test.set,x) ) 
  pred.test.nn=rbind(pred.test.nn, pred.nn )
  
  pred.bart=apply(as.matrix(bart_param),1,function(x) ML_list$bartfun(data=train.set,testdata=test.set, x) ) 
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

loss.function.gam=apply(pred.test.gam,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts.nn,tao))
loss.function.knn=apply(pred.test.knn,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts.nn,tao))
loss.function.svm=apply(pred.test.svm,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts.nn,tao))
loss.function.nn=apply(pred.test.nn,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts.nn,tao))
loss.function.rf=apply(pred.test.rf,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts.nn,tao))
loss.function.bart=apply(pred.test.bart,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts.nn,tao))
loss.function.lasso=apply(pred.test.lasso,2,function(x) ipcw_auc(T=data$ttilde,delta=data$delta,marker=crossprod(t(x),1),cause=1,wts=data$wts.nn,tao))

tuneparams=list(
  gam=grid.gam[which.max(loss.function.gam)],
  knn=grid.knn[which.max(loss.function.knn)],
  svm=grid.svm[which.max(loss.function.svm),],
  nn=grid.nn[which.max(loss.function.nn)],
  rf=grid.rf[which.max(loss.function.rf),],
  bart=grid.bart[which.max(loss.function.bart),],
  lasso=grid.lasso[which.max(loss.function.lasso)]
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

#svm: cost , gamma
#rf: number of trees, mtry
#bart: num_tree, k , q


save(tune.params,file="tune.params.RData")

