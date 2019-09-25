
# Simulations


### Functions ####
sim1 <- function(simdata,weighting,j){
  sim.data.train=simdata[[1]][[j]]
  sim.data.test=simdata[[2]][[j]]
  res <- ensBagg(train.data=sim.data.train,test.data=sim.data.test, xnam, tao=26.5 , weighting=weighting , folds=5,B=10 )
  true_ens <- apply(res$prediction_ensBagg,2, function(x) cvAUC::AUC(predictions =x,labels = sim.data.test$trueT))
  true_native <- apply(res$prediction_native_weights,2, function(x) cvAUC::AUC(predictions =x,labels = sim.data.test$trueT))
  true_survival <- apply(res$prediction_survival,2, function(x) cvAUC::AUC(predictions =x,labels = sim.data.test$trueT))
  return(rlist::list.append(res,true_ens=true_ens,true_native=true_native,true_survival=true_survival))
}

simulation1_all_scenarios <- function(weighting,d,s){
  simdata <- datagenPaper(J, n=1250 , frac.train=0.80 ,tao=26.5 , simulation=d, scenario=s )
  res_scen<- lapply(seq(1,J),function(x) sim1(simdata,weighting,x))
 }

table1_paper<- function(res,J,scenarios){
  num_scen<- length(scenarios)
  AUC_Boot_sim<- vector("list", num_scen)
  AUC_Native_sim <- vector("list", num_scen)
  AUC_Survival_sim <- vector("list", num_scen)
  AUC_true_ens_sim <- vector("list", num_scen)
  AUC_true_native_sim <- vector("list", num_scen)
  AUC_true_survival_sim <- vector("list", num_scen)
  table1 <- NULL
  for(s in 1:num_scen){
    for(j in 1:J){
      AUC_Boot_sim[[s]] <- rbind(AUC_Boot_sim[[s]],res[[s]][[j]]$auc_ipcwBagg)
      AUC_Native_sim[[s]] <- rbind(AUC_Native_sim[[s]],res[[s]][[j]]$auc_native_weights)
      AUC_Survival_sim[[s]] <- rbind(AUC_Survival_sim[[s]],res[[s]][[j]]$auc_survival)
      AUC_true_ens_sim[[s]] <- rbind(AUC_true_ens_sim[[s]],res[[s]][[j]]$true_ens)
      AUC_true_native_sim[[s]] <- rbind(AUC_true_native_sim[[s]],res[[s]][[j]]$true_native)
      AUC_true_survival_sim[[s]] <- rbind(AUC_true_survival_sim[[s]],res[[s]][[j]]$true_survival)
    }
    table1<- cbind(table1,
                   c(round(apply(AUC_Boot_sim[[s]],2,mean),3),
                     round(apply(AUC_Native_sim[[s]],2,mean),3),
                     round(apply(AUC_Survival_sim[[s]],2,mean),3)
                   ), c(round(apply(AUC_true_ens_sim[[s]],2,mean),3),
                        round(apply(AUC_true_native_sim[[s]],2,mean),3),
                        round(apply(AUC_true_survival_sim[[s]],2,mean),3)
                   )
    )
  }
  return(table1)
}

####     end of functions  #### 


J <- 500# number of simulations
xnam <- paste("X", 1:20, sep="") # names of the covariates 
scenarios=list(1,2,3,4)

#####    Weighting = CoxPH    ####

# Simulation 1
res_simulation1=lapply(scenarios, function(s) simulation1_all_scenarios(weighting = "CoxPH",d=1,s) )
table1_sim1 <- table1_paper(res_simulation1,J,scenarios)
# Simulation 2
res_simulation2=lapply(scenarios, function(s) simulation1_all_scenarios(weighting = "CoxPH",d=2,s) )
table1_sim2 <- table1_paper(res_simulation2,J,scenarios)

# Table 1: Average Estimated AUCs across 500 data sets for Simulation 1 and 2 and their four scenarios (A, B, C and D) using
# all available covariates and a Cox-PH model for censoring for predicting the event of interest in the test data set.

library(xtable)
xtable::xtable(cbind(table1_sim1,table1_sim2),digits=c(0,rep(3,16)))


#####    Weighting = CoxBoost   ####

# Simulation 1
res_simulation1=lapply(scenarios, function(s) simulation1_all_scenarios(weighting = "CoxBoost",d=1,s) )
table1_sim1 <- table1_paper(res_simulation1,J,scenarios)
# Simulation 2
res_simulation2=lapply(scenarios, function(s) simulation1_all_scenarios(weighting = "CoxBoost",d=2,s) )
table1_sim2 <- table1_paper(res_simulation2,J,scenarios)

# Table 2: Average Estimated AUCs across 500 data sets for Simulation 1 and 2 and their four scenarios (A, B, C and D) using
# all available covariates and a CoxBoost model for censoring for predicting the event of interest in the test data set.

xtable::xtable(cbind(table1_sim1,table1_sim2),digits=c(0,rep(3,16)))