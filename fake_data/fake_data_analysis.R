

################################################
#
# R script that runs the analysis of the fake data set. It is the same as used for the real data analysis. Analysis results cannot be exactly reproduced.
# We select the tuning parameters as it is described in the appendix of the paper  Gonzalez Ginestet et al. (2019+). "Ensemble IPCW Bagging bagging: a case study in the HIV care registry".
# The tuning parameter is done using the function ensBagg::tune_params_ml based on a grid values of each of the hyperparameter considered. 
# We consider the following hyperparameters:
# GAM: degree of freedom (df). we only consider df equal to 3 and 4
# LASSO: lambda
# Random Forest:  num_trees and mtry parameters.
# k-NN:  number of neighbours considered
# SVM: the cost parameter,  the gamma and the kernel. kernel=1 denotes "radial" and kernel=2 denotes "linear".
# Neural Network: number of neurons
# BartMachine:  num_tree parameter,  k parameter and  q parameter.
###############################################

rm(list = ls())
library(dplyr)
library(ensBagg)

#  we load the fake data set which lies in this folder
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


