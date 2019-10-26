

################################################
#
# R script that runs the analysis of the fake data set. It is the same as used for the real data analysis. Analysis results cannot be exactly reproduced.
# It reproduces:
# Table 3: Estimated AUCs for predicting a rebound in viral load in split sample test set based
# on Cox-PH-weights and Figure 1
# Figure 1: Aalen-Johansen estimates for the cause-specic cumulative
# incidences of failing to maintain HIV RNA undetectable within 2 years of suppression within
# quartiles of the test data set and  ROC curves for the the optimized ensemble IPCW bagged.
# 
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
                                                 data=fake_data_train,
                                                 weighting="CoxPH")

# we specify the algorithms
ens.library <-ensBagg::ens.all.algorithms()

# we predict the outcome of interest using the ensBagg with Cox PH weights using all algorithms available in the package
res.ensBagg <- ensBagg::ensBagg(train.data = fake_data_train,
                                test.data = fake_data_test,
                                xnam=xnam,
                                tao=730,
                                weighting = "CoxPH",
                                folds = 5,
                                ens.library = ens.library,
                                tuneparams = tuneparams_fakedata)


###########             Table 3                     ########### 
#Estimated AUCs for predicting a rebound in viral load in split sample test set based
#on Cox-PH-weights.

table3 <- cbind(c(res.ensBagg$auc_ipcwBagg,res.ensBagg$auc_survival[1:2]),c(res.ensBagg$auc_native_weights,rep(NA,7)),c(rep(NA,3),res.ensBagg$auc_survival[3],rep(NA,7)))
table3


######               Figure 1                       ##########
# Panel at right: the ROC curves for the the optimized ensemble IPCW bagged
ensBagg::plot_roc(time=fake_data_test$ttilde,delta=fake_data_test$delta, marker=res.ensBagg$prediction_ensBagg[,"Ensemble"], tao=730, method="ipcw")

#Panel at left: the Aalen-Johansen estimates for the cause-specific cumulative
#incidences of failing to maintain HIV RNA undetectable within 2 years of suppression within
#quartiles of the test data set

library(prodlim)
library(ggplot2)

tao=730

pred.ens <- res.ensBagg$prediction_ensBagg[,"Ensemble"]
fake_data_test$quartile <- with(fake_data_test,
                           cut(pred.ens, 
                               breaks=quantile(pred.ens, probs=seq(0,1, by=0.2), na.rm=TRUE), 
                               include.lowest=TRUE))
fake_data_test$quartile1 <- as.numeric(data.test$quartile)

cuminc_q=matrix(NA,max(fake_data_test$quartile1),1)
l_ci=matrix(NA,max(fake_data_test$quartile1),1)
u_ci=matrix(NA,max(fake_data_test$quartile1),1)
for(i in 1:max(fake_data_test$quartile1)){
 
  cum_inc_1 <- prodlim(Hist(ttilde,delta)~1,data=fake_data_test[fake_data_test$quartile1==i,],type = "cuminc")
  cuminc_q[i] <- cum_inc_1$cuminc$`1` [cum_inc_1$time>=tao][which.min(cum_inc_1$time>=tao)]
  
  u_ci[i]=cum_inc_1$upper$`1`[cum_inc_1$time>=tao][which.min(cum_inc_1$time>=tao)]
  l_ci[i]=cum_inc_1$lower$`1`[cum_inc_1$time>=tao][which.min(cum_inc_1$time>=tao)]
}

cuminc_all_test=prodlim(Hist(ttilde,delta)~1,data=fake_data_test,type = "cuminc")
ci_all.test=cuminc_all_test$cuminc$`1` [cuminc_all_test$time>=tao][which.min(cuminc_all_test$time>=tao)]
u.all.test=cuminc_all_test$upper$`1` [cuminc_all_test$time>=tao][which.min(cuminc_all_test$time>=tao)]
l.all.test=cuminc_all_test$lower$`1` [cuminc_all_test$time>=tao][which.min(cuminc_all_test$time>=tao)]

#data for the plot
df <- data.frame(x =levels(fake_data_test$quartile),
                 cuminc =cuminc_q,
                 L =l_ci,
                 U =u_ci)

df$x <- factor(df$x, levels = df$x[order(df$cuminc)])

figure1 <- ggplot(df, aes(x = factor(x), y = cuminc)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymax = U, ymin = L), width = 0.25)+
  scale_y_continuous(name="Cumulative Incidence at 2 years \n", breaks=c(0,.1,.2,.27,.3,.4,.5),limits=c(0, .5) )+
  xlab("\n Quartiles of IPCW bagging Super Learner Predictions")+
  theme_bw() + theme(panel.grid.major = element_blank(),plot.margin = margin(10,0,10,0),
                     axis.title.x = element_text( size=16),
                     axis.title.y = element_text( size=16),
                     axis.text.x = element_text( 
                       size=14),
                     axis.text.y = element_text( 
                       size=14)) +  theme(aspect.ratio = 1)+
  geom_hline(yintercept = ci_all.test) +
  geom_hline(yintercept = u.all.test,linetype="dotted")+
  geom_hline(yintercept = l.all.test,linetype="dotted")

figure1

