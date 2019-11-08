
# This script replicate the Simulation Study, section 5, from  Gonzalez Ginestet et al. (2019+). "Ensemble IPCW Bagging bagging: a case study in the HIV care registry"
# # It reproduces the figures showing the AUC of each single method relative to the ensemble and the tables with average estimated AUCs that are found in the supplementary material.
# The setting is the same as it is described in the paper:
# two ways of computing the IPC weights: Cox PH and Cox Boost
# 2 simulation studies and 4 scenarios in each simulation study
# n= 1250 sample size
# J= 500 simulated data sets
# 80% of the data is kept as training data set
# B=10 bootstrap samples
# folds=5 number of folds
# tune parameters are those that are specified in the paper and were computed using the function EnsBagg::parametersSimulation(folds = 5,xnam,train.data,tao) that computes the default tune parameters.
# Expected run time (standard desktop) of one simulation study (which is 4 scenarios and 500 simulated data sets) is 235 hours (9.8 days)  


### Functions ####
sim <- function(simdata,weighting,j){
  sim.data.train=simdata[[1]][[j]][,-(1:2)] # remove the first two columns: id and E so to have the appropiate format
  sim.data.test=simdata[[2]][[j]][,-(1:2)] # remove the first two columns: id and E
  res <- ensBagg(train.data=sim.data.train,test.data=sim.data.test, xnam, tao=26.5 , weighting=weighting , folds=5,ens.library=ens.library ,B=10 )
  true_ens <- apply(res$prediction_ensBagg,2, function(x) cvAUC::AUC(predictions =x,labels = sim.data.test$trueT))
  true_native <- apply(res$prediction_native_weights,2, function(x) cvAUC::AUC(predictions =x,labels = sim.data.test$trueT))
  true_survival <- apply(res$prediction_survival,2, function(x) cvAUC::AUC(predictions =x,labels = sim.data.test$trueT))
  return(rlist::list.append(res,true_ens=true_ens,true_native=true_native,true_survival=true_survival,auc_true=simdata[[3]][[j]]))
}

simulation_all_scenarios <- function(weighting,d,s){
  simdata <- datagenPaper(J, n=1250 , frac.train=0.80 ,tao=26.5 , simulation=d, scenario=s )
  res_scen<- lapply(seq(1,J),function(x) sim(simdata,weighting,x))
 }

table1_paper<- function(res,J,scenarios){
  num_scen<- length(scenarios)
  AUC_Boot_sim<- vector("list", num_scen)
  AUC_Native_sim <- vector("list", num_scen)
  AUC_Survival_sim <- vector("list", num_scen)
  AUC_true_ens_sim <- vector("list", num_scen)
  AUC_true_native_sim <- vector("list", num_scen)
  AUC_true_survival_sim <- vector("list", num_scen)
  AUC_true <- vector("list", num_scen)
  table1 <- NULL
  for(s in 1:num_scen){
    for(j in 1:J){
      AUC_Boot_sim[[s]] <- rbind(AUC_Boot_sim[[s]],res[[s]][[j]]$auc_ipcwBagg)
      AUC_Native_sim[[s]] <- rbind(AUC_Native_sim[[s]],res[[s]][[j]]$auc_native_weights)
      AUC_Survival_sim[[s]] <- rbind(AUC_Survival_sim[[s]],res[[s]][[j]]$auc_survival)
      AUC_true_ens_sim[[s]] <- rbind(AUC_true_ens_sim[[s]],res[[s]][[j]]$true_ens)
      AUC_true_native_sim[[s]] <- rbind(AUC_true_native_sim[[s]],res[[s]][[j]]$true_native)
      AUC_true_survival_sim[[s]] <- rbind(AUC_true_survival_sim[[s]],res[[s]][[j]]$true_survival)
      AUC_true[[s]] <- rbind(AUC_true[[s]],res[[s]][[j]]$auc_true)
    }
    table1<- cbind(table1,
                   c(round(apply(AUC_true[[s]],2,mean),3),
                     round(apply(AUC_Boot_sim[[s]],2,mean),3),
                     round(apply(AUC_Native_sim[[s]],2,mean),3),
                     round(apply(AUC_Survival_sim[[s]],2,mean),3)
                   ), c(NA,
                     round(apply(AUC_true_ens_sim[[s]],2,mean),3),
                        round(apply(AUC_true_native_sim[[s]],2,mean),3),
                        round(apply(AUC_true_survival_sim[[s]],2,mean),3)
                   )
    )
  }
  return(table1)
}


figure1_paper<- function(res,J,scenarios){
  num_scen<- length(scenarios)
  AUC_sim<- vector("list", num_scen)
  data_boxplot<- vector("list", num_scen)
  upper_panel <- matrix(NA,nrow = J,17)
  for(j in 1:J){
    upper_panel[j,] <- c(rep(1,9),rep(2,5),rep(3,3))
  }
  colnames(upper_panel) <- c(rep("IPCW Bagging",9),rep("Native Weight",5),rep("Survival",3))
 
  for(s in 1:num_scen){
    for(j in 1:J){
      
      AUC_sim[[s]] <- rbind(AUC_sim[[s]],c(res[[s]][[j]]$auc_ipcwBagg[,-10],res[[s]][[j]]$auc_native_weights,res[[s]][[j]]$auc_survival)/res[[s]][[j]]$auc_ipcwBagg[,10] )
      
    }
    colnames(AUC_sim[[s]]) <- c("LR","GAM.3","GAM.4","LASSO","RF","SVM","BART","k-NN","NNet","LR","GAM.3","GAM.4",
                                "LASSO",  "RF", "Cox" ,"CoxBoost" , "RF")
    data_boxplot[[s]] <- melt(AUC_sim[[s]])
    data_boxplot[[s]]$scenario <- 1
    data_boxplot[[s]] <- cbind(data_boxplot[[s]], melt(upper_panel))
    
  }
  
  data_boxplot1<- rbind( data_boxplot[[1]],data_boxplot[[2]],data_boxplot[[3]],data_boxplot[[4]] )
  data_boxplot1 <- data_boxplot1[,c(-1,-5,-7)]
  colnames(data_boxplot1) <- c("ml","r.auc","scenario","method")
  
  p <- ggplot(data_boxplot1, aes(ml, r.auc,fill=factor(scenario)))
  p <- p + geom_boxplot()  + labs(x="\n Machine Learning (ML) Algorithms",y=TeX("$AUC_{ML} /  AUC_{Ensemble} \n "),fill = "Scenario") +
    facet_grid(~method,  scales = "free_x", space = "free_x") +
    geom_hline(aes(yintercept = 1)) + scale_y_continuous(breaks=seq(0.5,1.2,.05)) +
    theme_bw() +  theme(panel.grid.major = element_blank(),plot.margin = margin(10,0,10,0),
                        axis.title.x = element_text( size=25),
                        axis.title.y = element_text( size=25),
                        axis.text.x = element_text( 
                          size=20),
                        axis.text.y = element_text( 
                          size=20)) +
    theme(strip.text.x = element_text(size = 25)) +
    theme(legend.position = "none") +  scale_fill_grey()
  
  return(p)
  
}

####     end of functions  #### 
#libraries for the plot
library(reshape2)
library(ggplot2)
library(tikzDevice)
library(latex2exp)

J <- 500# number of simulations
xnam <- paste("X", 1:20, sep="") # names of the covariates 
scenarios=list(1,2,3,4)
ens.library <- ensBagg::ens.all.algorithms()
#####    Weighting = CoxPH    ####

# Simulation 1
res_simulation1=lapply(scenarios, function(s) simulation_all_scenarios(weighting = "CoxPH",d=1,s) )
table1_sim1 <- table1_paper(res_simulation1,J,scenarios)
row.names(table1_sim1)[1] <- "True"
figure1_sim1 <-figure1_paper(res_simulation1,J,scenarios) 

# Simulation 2
res_simulation2=lapply(scenarios, function(s) simulation_all_scenarios(weighting = "CoxPH",d=2,s) )
table1_sim2 <- table1_paper(res_simulation2,J,scenarios)
row.names(table1_sim2)[1] <- "True"
figure1_sim2 <-figure1_paper(res_simulation2,J,scenarios) 


# AUCs of each algorithm relative to the AUC of the ensemble based on 500 simulated data sets
#under the four scenarios (A being the darkest and D the whitest, alphabetically order)  using
#all available covariates and a Cox-PH model for censoring for predicting the event of interest in the test data set.
#The horizontal line denotes the ensemble

figure1_sim1 #simulation 1
figure1_sim2 #simulation 2


# Table 1: Average Estimated AUCs across 500 data sets for Simulation 1 and 2 and their four scenarios (A, B, C and D) using
# all available covariates and a Cox-PH model for censoring for predicting the event of interest in the test data set.

library(xtable)
xtable::xtable(cbind(table1_sim1,table1_sim2),digits=c(0,rep(3,16)))


#####    Weighting = CoxBoost   ####

# Simulation 1
res_simulation1=lapply(scenarios, function(s) simulation1_all_scenarios(weighting = "CoxBoost",d=1,s) )
table1_sim1 <- table1_paper(res_simulation1,J,scenarios)
row.names(table1_sim1)[1] <- "True"
figure1_sim1 <-figure1_paper(res_simulation1,J,scenarios) 

# Simulation 2
res_simulation2=lapply(scenarios, function(s) simulation1_all_scenarios(weighting = "CoxBoost",d=2,s) )
table1_sim2 <- table1_paper(res_simulation2,J,scenarios)
row.names(table1_sim2)[1] <- "True"
figure1_sim2 <-figure1_paper(res_simulation1,J,scenarios) 

# AUCs of each algorithm relative to the AUC of the ensemble based on 500 simulated data sets
#under the four scenarios (A being the darkest and D the whitest, alphabetically order)  using
#all available covariates and a CoxBoost model for censoring for predicting the event of interest in the test data set.
#The horizontal line denotes the ensemble

figure1_sim1 #simulation 1
figure1_sim2 #simulation 2


# Table 2: Average Estimated AUCs across 500 data sets for Simulation 1 and 2 and their four scenarios (A, B, C and D) using
# all available covariates and a CoxBoost model for censoring for predicting the event of interest in the test data set.

xtable::xtable(cbind(table1_sim1,table1_sim2),digits=c(0,rep(3,16)))