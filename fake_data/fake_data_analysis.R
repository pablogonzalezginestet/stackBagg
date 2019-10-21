


rm(list = ls())
library(dplyr)
library(ensBagg)

#  we load the fake data set
load("fake_data_train.RData")
load("fake_data_test.RData")

# we subset the train fake data set because the random forest crashes  but one could run all the following 
# lines using the 2335 observations of fake_data_train 

fake_data_train <- fake_data_train[sample(nrow(fake_data_train),1000,replace = FALSE),]

# we define the variables to be included in the analysis
xnam <- names(fake_data_train)[-(1:2)]

# we set the grid of the hyperparameters that we are going to use in the tuning step
grid.hyperparam <- ensBagg::grid_parametersDataHIV(xnam,fake_data_train,tao=730)

# we tune the hyperparameters using the grid 
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

# we specify the algorithms
ens.library <-ensBagg::ens.all.algorithms()

# we predict the outcome of interest using the ensBagg with Cox PH weights using all algorithms available in the package
res.ensBagg <- ensBagg::ensBagg(train.data = fake_data_train,test.data = fake_data_test,xnam=xnam,tao=730,weighting = "CoxPH",folds = 5,ens.library = ens.library )


