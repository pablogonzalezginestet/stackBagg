
# Simulations


########### Functions #############
sim1 <- function(simdata, j){
  sim.data.train=simdata[[1]][[j]]
  sim.data.test=simdata[[2]][[j]]
  res <- ensBagg(train.data=sim.data.train,test.data=sim.data.test, xnam, tao=26.5 , weighting="CoxPH" , folds=2,B=2 )
  true_ens <- apply(res$prediction_ensBagg,2, function(x) cvAUC::AUC(predictions =x,labels = sim.data.test$trueT))
  true_native <- apply(res$prediction_native_weights,2, function(x) cvAUC::AUC(predictions =x,labels = sim.data.test$trueT))
  true_survival <- apply(res$prediction_survival,2, function(x) cvAUC::AUC(predictions =x,labels = sim.data.test$trueT))
  return(rlist::list.append(res,true_ens=true_ens,true_native=true_native,true_survival=true_survival))
}

simulation1_all_scenarios <- function(d,s){
  simdata <- datagenPaper(J, n=500 , frac.train=0.80 ,tao=26.5 , simulation=d, scenario=s )
  res_scen<- lapply(seq(1,J),function(x) sim1(simdata,x))
 }

##################################

J <- 500 # number of simulations
xnam <- paste("X", 1:20, sep="") # names of the covariates 
scenarios=list(1,2,3,4)
res_simulation1=lapply(scenarios, function(s) simulation1_all_scenarios(d=1,s) )

#########################################################################


AUC_Boot_sim1<- vector("list", 4)
AUC_Native_sim1 <- vector("list", 4)
AUC_Survival_sim1 <- vector("list", 4)
AUC_true_ens_sim1 <- vector("list", 4)
AUC_true_native_sim1 <- vector("list", 4)
AUC_true_survival_sim1 <- vector("list", 4)
for(s in 1:4){
for(j in 1:J){
  AUC_Boot_sim1[[s]] <- rbind(AUC_Boot_sim1[[s]],res_simulation1[[s]][[j]]$auc_ipcwBagg)
  AUC_Native_sim1[[s]] <- rbind(AUC_Native_sim1[[s]],res_simulation1[[s]][[j]]$auc_native_weights)
  AUC_Survival_sim1[[s]] <- rbind(AUC_Survival_sim1[[s]],res_simulation1[[s]][[j]]$auc_survival)
  AUC_true_ens_sim1[[s]] <- rbind(AUC_true_ens_sim1[[s]],res_simulation1[[s]][[j]]$true_ens)
  AUC_true_native_sim1[[s]] <- rbind(AUC_true_native_sim1[[s]],res_simulation1[[s]][[j]]$true_native)
  AUC_true_survival_sim1[[s]] <- rbind(AUC_true_survival_sim1[[s]],res_simulation1[[s]][[j]]$true_survival)
    }
}



table1_cox_weights <- NULL
for(s in 1:4){
  table1_cox_weights<- cbind(table1_cox_weights,
    c(
    round(apply(AUC_Boot_sim1[[s]],2,mean),3),
    round(apply(AUC_Native_sim1[[s]],2,mean),3),
    round(apply(AUC_Survival_sim1[[s]],2,mean),3)
  ),
   c(round(apply(AUC_true_ens_sim1[[s]],2,mean),3),
    round(apply(AUC_true_native_sim1[[s]],2,mean),3),
    round(apply(AUC_true_survival_sim1[[s]],2,mean),3)
  )
  )
}


library(xtable)
