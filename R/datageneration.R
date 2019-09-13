
# J = number of simulated data sets
# n = sample size
# frac.train = percentange train data set. A number between 0 and 1.
# simulation = simulation indicator setting. It takes on 1 and 2. 
# scenario = scenario indicator. It takes on 1, 2, 3 and 4.

# scenario=1 : Independent censoring. Censoring time does not depend on any covariate.
# scenario=2 : Independent censoring, non-informative. Censoring time depends on a
# disjoint subset of the covariates from those associated with the failure outcome.
# Scenario=3 : Informative, dependent censoring. Censoring time depends on the same 
#subset of the covariates from those associated with the failure outcome
# scenario=4 : Informative, dependent censoring. Censoring time and the competing
#event time depend on the same subset of the covariates from those associated with the
#failure outcome. 

# g1 = scale parameter type event 1
# g2 = shape parameter type event 1
# k1 = scale parameter type event 2
# k2 = shape parameter type event 2

# cenper = percentage level of censoring
# typ1expt = percertange of the individuals suffering type 1 event
# typ2expt = percertange of the individuals suffering type 1 event

datagenPaper=function(J = 500 , n = 1250, frac.train = 0.8, simulation, scenario) {
  
  set.seed(500)
  train.data=list()
  test.data=list()
  Tstar=26.5
  tao= Tstar
  D <- matrix(NA,nrow = J,ncol = 3)
  
  for (j in 1:J){
    
    #############  Data Generating Process ##############
    X5 <- matrix(rnorm(n * 5, mean = 0, sd = .1), ncol = 5)
    X52 <- X5 %*% matrix(.25, nrow = 5, ncol = 5) + matrix(rnorm(n * 5, mean = 0, sd = .1), ncol = 5)
    X53 <- X52 %*% matrix(.15, nrow = 5, ncol = 5) + matrix(rnorm(n * 5, mean = 0, sd = .5), ncol = 5)
    X54 <- X53 %*% matrix(.05, nrow = 5, ncol = 5) + matrix(rnorm(n * 5, mean = 0, sd = .65), ncol = 5)
    
    X <- cbind(X5, X52, X53, X54)
    
    
    if(simulation==1){

    X2 <- X[, c(1, 6, 11, 16)]
    X2 <- cbind(X2[,1],X2[,2], X2[, 1] * X2[, 2], bs(X2[, 3], df = 4, degree = 3), bs(X[, 4], df = 4, degree = 3))
    
    beta <- c(.75, .75, 5, 1.5, -2.5, 3, 2, 2, 3, 3, .9)
    g1 <- sqrt(exp(6 + X2 %*% beta )) 
    g2 <- 3.5  
    
    if(scenario == 1) {
      ###############   SCENARIO A  #################
      
      k1 <- exp(2.5 * X[, 2]) 
      k2 <- 2.5 
      cenper<-0.22 
      cskl <- (Tstar / ((-log(1 - cenper)) ^ (1)))
      
      typ1expt<- 0.27  
      typ2expt<- 0.1 
      rescl <- (Tstar / ((-log(1-typ1expt)) ^ (1/g2)))
      reskl <- (Tstar / ((-log(1-typ2expt)) ^ (1 / k2)))
      
      g1 <- g1 * mean(rescl / (g1))
      k1 <- k1 * mean(reskl / (k1))
      
      Y <- rweibull(n, scale = g1, shape = g2)
      Y2 <- rweibull(n, scale = k1, shape = k2)
      Cen <- rweibull(n, scale = cskl, shape = 1)
      
      
    } else if(scenario == 2) {
      ###############   SCENARIO B  #################
      
      k1 <- exp(2.5 * X[, 2]) 
      k2 <- 2.5  
      
      censb1 <- sqrt(exp( X[, c(3, 4,9 ,19, 20)] %*% c(.7,.7,.7,-2,-2)))
      cenper<-0.26
      cskl <- (Tstar / ((-log(1 - cenper)) ^ (1)))
      censb1 <- censb1 * mean(cskl / censb1)
      
      typ1expt<- 0.27
      typ2expt<- 0.1 
      
      rescl <- (Tstar / ((-log(1-typ1expt)) ^ (1/g2)))
      reskl <- (Tstar / ((-log(1-typ2expt)) ^ (1 / k2)))
      
      g1 <- g1 * mean(rescl / (g1))
      k1 <- k1 * mean(reskl / (k1))
      
      Y <- rweibull(n, scale = g1, shape = g2)
      Y2 <- rweibull(n, scale = k1, shape = k2)
      Cen <- rweibull(n, scale = censb1, shape = 1)
      
    } else if(scenario == 3) {
      
      k1 <- exp(2.5 * X[, 2]) 
      k2 <- 2.5
      
      censb1 <- sqrt(exp( X[, c(1, 6, 11, 16,20)] %*% c(.7,.7,.7,-2,-2)))
      cenper<-0.26
      cskl <- (Tstar / ((-log(1 - cenper)) ^ (1)))
      censb1 <- censb1 * mean(cskl / censb1)
      
      typ1expt<- 0.27  
      typ2expt<- 0.1 
      
      rescl <- (Tstar / ((-log(1-typ1expt)) ^ (1/g2)))
      reskl <- (Tstar / ((-log(1-typ2expt)) ^ (1 / k2)))
      
      g1 <- g1 * mean(rescl / (g1))
      k1 <- k1 * mean(reskl / (k1))
      
      Y <- rweibull(n, scale = g1, shape = g2)
      Y2 <- rweibull(n, scale = k1, shape = k2)
      Cen <- rweibull(n, scale = censb1, shape = 1)
      
    } else if(scenario == 4) {
      
      k1 <- exp( X[, c(1, 6, 11, 16,20)] %*% c(-1,-.9,1,2,2))
      k2 <- 2.5 
      
      censb1 <- sqrt(exp( X[, c(1, 6, 11, 16,20)] %*% c(.7,.7,.7,-2,-2)))
      cenper<-0.26
      cskl <- (Tstar / ((-log(1 - cenper)) ^ (1)))
      censb1 <- censb1 * mean(cskl / censb1)
 
      typ1expt<- 0.27  
      typ2expt<- 0.1 
      
      rescl <- (Tstar / ((-log(1-typ1expt)) ^ (1/g2)))
      reskl <- (Tstar / ((-log(1-typ2expt)) ^ (1 / k2)))
      
      g1 <- g1 * mean(rescl / (g1))
      k1 <- k1 * mean(reskl / (k1))
      
      Y <- rweibull(n, scale = g1, shape = g2)
      Y2 <- rweibull(n, scale = k1, shape = k2)
      Cen <- rweibull(n, scale = censb1, shape = 1)
      
    }
    ######  simulation==2
    }else{
      
      X2 <- X[, c(1, 6, 11, 16, 20)]
      X2 <- cbind(X2, X2[, 1] * X2[, 2],
                  cos(X2[, 3] / .1),
                  X2[, 4] * ifelse(X2[, 4] < median(X2[, 4]), 0, 1))
      
      
      beta.c <- c(1.1, 1.4, -2.1, -1.2, -2.3, -1.5, 6.7, .5) / 4
      g1 <- sqrt(exp(X2 %*% beta.c))
      g2 <- 3.5 
      
      if(scenario == 1) {
        ###############   SCENARIO A  #################
        
        k1 <- exp(2.5 * X[, 2]) 
        k2 <- 2.5 
        cenper<-0.22 
        cskl <- (Tstar / ((-log(1 - cenper)) ^ (1)))
        
        typ1expt<- 0.24 
        typ2expt<- 0.1  
        
        rescl <- (Tstar / ((-log(1-typ1expt)) ^ (1/g2)))
        reskl <- (Tstar / ((-log(1-typ2expt)) ^ (1 / k2)))
        
        g1 <- g1 * mean(rescl / (g1))
        k1 <- k1 * mean(reskl / (k1))
        
        Y <- rweibull(n, scale = g1, shape = g2)
        Y2 <- rweibull(n, scale = k1, shape = k2)
        Cen <- rweibull(n, scale = cskl, shape = 1)
        
        
      } else if(scenario == 2) {
        ###############   SCENARIO B  #################
        
        k1 <- exp(2.5 * X[, 2]) 
        k2 <- 2.5 
       
        censb1 <- sqrt(exp( cbind(ifelse(X[,3]<median(X[,3]),0,1),X[,4], X[,9],cos(X[,19])/.7) %*% c(.5,.7,.7,.8)))
        cenper<-0.26
        cskl <- (Tstar / ((-log(1 - cenper)) ^ (1)))
        censb1 <- censb1 * mean(cskl / censb1)
      
        
        typ1expt<- 0.24  
        typ2expt<- 0.1 
        
        rescl <- (Tstar / ((-log(1-typ1expt)) ^ (1/g2)))
        reskl <- (Tstar / ((-log(1-typ2expt)) ^ (1 / k2)))
        
        g1 <- g1 * mean(rescl / (g1))
        k1 <- k1 * mean(reskl / (k1))
        
        Y <- rweibull(n, scale = g1, shape = g2)
        Y2 <- rweibull(n, scale = k1, shape = k2)
        Cen <- rweibull(n, scale = censb1, shape = 1)
        
      } else if(scenario == 3) {
        
        k1 <- exp(2.5 * X[, 2])
        k2 <- 2.5  
        
        censb1 <- sqrt(exp( cbind(ifelse(X[,1]<median(X[,1]),0,1),X[,6], X[,11],cos(X[,16])/.7,X[,20]) %*% c(.5,.7,.7,.8,-2)))
        cenper<-0.26
        cskl <- (Tstar / ((-log(1 - cenper)) ^ (1)))
        censb1 <- censb1 * mean(cskl / censb1)
        
        
        typ1expt<- 0.24
        typ2expt<- 0.1 
     
        rescl <- (Tstar / ((-log(1-typ1expt)) ^ (1/g2)))
        reskl <- (Tstar / ((-log(1-typ2expt)) ^ (1 / k2)))
        
        g1 <- g1 * mean(rescl / (g1))
        k1 <- k1 * mean(reskl / (k1))
        
        Y <- rweibull(n, scale = g1, shape = g2)
        Y2 <- rweibull(n, scale = k1, shape = k2)
        Cen <- rweibull(n, scale = censb1, shape = 1)
        
        
      } else if(scenario == 4) {
        
        k1 <- exp( X[, c(1, 6, 11, 16,20)] %*% c(-1,-.9,1,2,2))
        k2 <- 2.5 
        
        censb1 <- sqrt(exp( cbind(ifelse(X[,1]<median(X[,1]),0,1),X[,6], X[,11],cos(X[,16])/.7,X[,20]) %*% c(.5,.7,.7,.8,-2)))
        cenper<-0.26
        cskl <- (Tstar / ((-log(1 - cenper)) ^ (1)))
        censb1 <- censb1 * mean(cskl / censb1)
        
        typ1expt<- 0.24 
        typ2expt<- 0.1 
        
        rescl <- (Tstar / ((-log(1-typ1expt)) ^ (1/g2)))
        reskl <- (Tstar / ((-log(1-typ2expt)) ^ (1 / k2)))
        
        g1 <- g1 * mean(rescl / (g1))
        k1 <- k1 * mean(reskl / (k1))
        
        Y <- rweibull(n, scale = g1, shape = g2)
        Y2 <- rweibull(n, scale = k1, shape = k2)
        Cen <- rweibull(n, scale = censb1, shape = 1)
        
      }
      
      
    }
    
    ttilde <- pmin(Y, Y2, Cen)
    delta <- ifelse(Cen < Y & Cen < Y2, 0,
                    ifelse(Y < Y2, 1, 2))
    
    
    sim.data=data.frame(id = 1:length(ttilde), ttilde, delta, trueT = Y < Tstar & Y < Y2, Cen, Y, Y2,X)
    sim.data<- mutate(sim.data,E=as.factor(ifelse(ttilde < Tstar & delta==1, 1 , ifelse(ttilde < Tstar & delta==2 | ttilde>Tstar, 0, NA))),
                      deltac.1=ifelse(delta==0,1,0))
    
    xnam <- paste("X", 1:20, sep="")

    train.set <- sample(1:nrow(sim.data), floor(frac.train*nrow(sim.data)), replace=FALSE) 
    test.set <- setdiff(1:nrow(sim.data), train.set) 
    sim.data.train <- data.frame(sim.data[train.set,])
    sim.data.test <- data.frame(sim.data[test.set,]) 
    
    sim.data.train<- sim.data.train[c("id","E","wts.1","wts.11","ttilde","delta","trueT",xnam)]
    colnames(sim.data.train)[c(3,4)] <- c("wts","wts.nn")
    sim.data.test.1 <- sim.data.test[c("id","E","wts.1","wts.11","ttilde","delta","trueT",xnam)]
    colnames(sim.data.test.1)[c(3,4)] <- c('wts',"wts.nn")
    
    train.data[[j]]=sim.data.train
    test.data[[j]]=sim.data.test
    
  }  
  
  return(list(train.data, test.data))
  
}
