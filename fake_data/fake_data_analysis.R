


rm(list = ls())
library(dplyr)
library(ensBagg)

load("fake_data_train.RData")
load("fake_data_test.RData")


head(train)
head(test)

xnam <- names(fake_data_train)[-(1:2)]

fake_data_train <- fake_data_train[sample(nrow(fake_data_train),1000,replace = FALSE),]

grid.hyperparam <- ensBagg::grid_parametersDataHIV(xnam,fake_data_train,tao=730)

tuneparams_fakedata<- ensBagg::tune_params_ml(gam_param = grid.hyperparam$gam_param, 
                                                 lasso_param = grid.hyperparam$lasso_param,
                                                 randomforest_param = grid.hyperparam$randomforest_param,
                                                 knn_param = grid.hyperparam$knn_param,
                                                 svm_param = grid.hyperparam$svm_param,
                                                 nn_param = grid.hyperparam$nn_param,
                                                 bart_param = grid.hyperparam$bart_param,
                                                 folds=5,
                                                 xnam=xnam,
                                                 tao=730,
                                                 data=fake_data_train )

pred_fakedata<- ensBagg::ensBagg(train.data = train,test.data = test,xnam = xnam,tao = 730,weighting = "CoxBoost",folds = 5,tuneparams=tuneparams_fakedata)



