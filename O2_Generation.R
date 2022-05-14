library(survival)
library(copula)
library(doParallel)
library(randomizr)
library(ggplot2)
registerDoParallel(detectCores())
set.seed(3)




## Parameters for O2 parallel cluster
## inputTau = correlation between outcomes within each cluster
## inputClusterType = prespecified cluster types
# strataStatus = turn on stratified data generation
# randomization = type of randomization
params <- commandArgs(trailingOnly=TRUE)
inputTau <- as.numeric(params[[1]]) #0, 0.1, 0.2, 0.3, 0.7
inputClusterType <- as.numeric(params[[2]]) #1=small, 2=large, 3=mixed
strataStatus <- as.numeric(params[[3]]) #0=FALSE, 1=TRUE
randomization <- as.character(params[[4]]) #Individual, Cluster, Within-cluster

## Parameters for data generation
# simCount = number of simulations
# n = total number of people
# censor.rate = exponential rate for censor time generation
# beta = parameter of interest
simCount<-3000
n<-30000
censor.rate<-0.225
beta<-log(0.4) 




## Main data generator
# Generates correlated outcome data using a copula clayton model - Corresponds to a gamma frailty model with shape/scale=1
# n = total number of people
# K = cluster size 
# tau = correlation between outcomes within each cluster
# beta = parameter of interest
# randomization = type of randomization: "Individual", "Cluster", "Within-cluster" 
# strata = turn on stratified data generation: T/F
# censor.rate = exponential rate for censor time generation
# custom_nk = custom nk vector for mixed cluster sizes
generateData.copula <- function(n, K, tau, beta, randomization, strata, censor.rate, custom_nk=NULL){
  # Number of people per cluster 
  if (!is.null(custom_nk)) # Pass in a custom nk vector
    nk <- custom_nk
  else # Equally sized clusters
    nk <- rep(K, n/K)
  
  # Randomly generate covariate vector X for each person
  if (randomization=="Individual") # Individual randomization
    X <- complete_ra(N = n, m = n/2)
  else if (randomization=="Cluster") # Cluster randomization
    X <- cluster_ra(clusters=rep(1:length(nk), nk))
  else if (randomization=="Within-cluster") # Blocked randomization
    X <- block_ra(blocks=rep(1:length(nk), nk))
  
  # Copula parameter
  cpl.param <- 2*tau/(1-tau)
  
  # Initialize data frame
  mydat <- data.frame(matrix(ncol = 5, nrow = n))
  colnames(mydat) <- c("grp", "grp_size", "X", "U", "delta")
  curr_cluster <- 1
  start <- 1 # index for first person in current cluster
  end <- 0 # index for last person in current cluster
  
  # Generate the data
  for (i in nk){
    end <- end + i
    
    # Generate S - marginal survival functions s(T|X)
    if (tau!=0){ # correlated outcomes for tau > 0
      # Construct a copula object with alpha=cpl.param and dimension=current cluster size (nk[i])
      cpl_type <- claytonCopula(param = cpl.param, dim = i)
      S <- as.numeric(rCopula(1,  cpl_type))
    }
    else{
      S <- runif(i) # independent outcomes when tau = 0
    }
    
    if (strata==TRUE){
      # Method 1
      # For each cluster, set a different baseline hazard by uniformly selecting shape
      stratify <- sample(1:10, 1)
      t <- qweibull(exp(log(S)/exp(beta*X[start:end])), shape = stratify, scale = 70, lower.tail = F)
    }
    else{
      # Generate T - time to event
      t <- qweibull(exp(log(S)/exp(beta*X[start:end])), shape = 2, scale = 70, lower.tail = F )
    }
    
    # Generate C - censor time
    C <- rexp(n = i, rate = censor.rate)
    
    # Generate delta - indicator for whether a person is censored or not
    mydat[start:end,"delta"] <- ifelse(t <= C, 1, 0)
    
    # Save the minimum of the two times as the recorded time
    mydat[start:end,"U"] <- pmin(t, C)
    
    # Record cluster group and cluster size
    mydat[start:end,"grp"] <- rep(curr_cluster, i)
    mydat[start:end,"grp_size"] <- rep(i , i)
    
    # Increment indices
    curr_cluster <- curr_cluster+1
    start <- end+1
  }
  
  # Add covariate X to data frame
  mydat[,"X"] <- X
  
  return(mydat) 
}



## Helper Function
# Generates naive and robust confidence intervals, given a data set
genCI <- function(input_data, strata=FALSE){
  # Fit a cox proportional hazards model with robust option on
  # U = Time to event outcome, delta = censor variable, X = covariate of interest, grp = cluster group number
  if (strata==TRUE)
    fit <- coxph(Surv(U, delta) ~ X + strata(grp), cluster = grp, data = input_data)
  else
    fit <- coxph(Surv(U, delta) ~ X, cluster = grp, data = input_data)
  
  # Obtain se, robust se, and beta_hat
  beta.se <- sqrt(fit$naive.var[1])
  beta.rse <- sqrt(fit$var[1])
  beta.est <-fit$coefficients
  estimates <- c(est=beta.est, se=beta.se, rse=beta.rse)
  
  # Obtain the confidence intervals
  robust_ci <- confint(fit)
  naive_ci <- fit$coefficients[1] + c(-1,1)*1.96*beta.se

  all_ci<-list(estimates=estimates, naive=as.numeric(naive_ci), robust=as.numeric(robust_ci))
  return(all_ci)
}



## Helper Function
# Calculate width and coverage of CIs, given the output of runSimulation and true value of the parameter
calcCIProperties <- function(input_ci){
  # Calculate the widths of the different CIs
  width <- lapply(input_ci, function(x){
    naive_width <- abs(x[["naive"]][1] - x[["naive"]][2])
    robust_width <- abs(x[["robust"]][1] - x[["robust"]][2])
    c(naive_width=naive_width, robust_width=robust_width)
  })
  
  # Convert to data frame
  width <- do.call(rbind, width)
  average_width <- apply(width, 2, mean)
  
  # Calculate coverage
  naive_sum<-robust_sum<-rand_sum<-0
  for (i in input_ci){
    if (beta >= i[["naive"]][1] && beta <= i[["naive"]][2])
      naive_sum<-naive_sum+1
    if (beta >= i[["robust"]][1] && beta <= i[["robust"]][2])
      robust_sum<-robust_sum+1
  }
  
  simCount<-length(input_ci)
  coverage<-c(naive=naive_sum/simCount, robust=robust_sum/simCount)
  
  return(list(avg_width=average_width, coverage=coverage))
}




## Helper Function
# Calculates and returns data frame of statistics of interest
calcStats <- function(results, strata, genCI){
  # Calculate beta_hat, naive and robust s.e., and confidence intervals for each simulation
  ci.copula <- foreach(i=1:simCount, .packages = 'survival') %dopar% genCI(results[[i]], strata)
  
  # Eliminate very extreme values (>5 sd away)
  beta.hat.values <- sapply(ci.copula, function(x){x$estimates[[1]]})
  ese.robust <- mean(sapply(ci.copula, function(x){x$estimates[[3]]}))
  extreme.indices <- which(abs(beta.hat.values-beta) > (ese.robust*5))
  num.removed <- length(extreme.indices)
  
  if (num.removed > 0)
    ci.copula <- ci.copula[ (-extreme.indices) ]

  # Calculate beta_hat, empirical SE, and average estimated SE for naive and robust methods
  beta.hat.values <- sapply(ci.copula, function(x){x$estimates[[1]]})
  beta.hat <- mean(beta.hat.values)
  se.empirical <- sd(sapply(ci.copula, function(x){x$estimates[[1]]}))
  ese.naive <- mean(sapply(ci.copula, function(x){x$estimates[[2]]}))
  event.mean <- mean(sapply(results.copula, function(x){sum(x$delta)}))

  ci.properties<-calcCIProperties(ci.copula)
  
  # Plot histograms of beta.hat values
  plot <- qplot(beta.hat.values, geom="histogram", binwidth=0.05)
  if (strata == FALSE)
    ggsave(plot,file=paste0("Tau", inputTau, "Cluster", clusterString, "Random", randomization, "Strata", strataStatus,"_unstratifiedStats", ".pdf") )
  else
    ggsave(plot,file=paste0("Tau", inputTau, "Cluster", clusterString, "Random", randomization, "Strata", strataStatus, "_stratifiedStats", ".pdf") )
  
  final.table <- data.frame(beta.hat=beta.hat, se.empirical=se.empirical, ese.naive=ese.naive, cp.naive=ci.properties[[2]][1],
                            width.naive=ci.properties[[1]][1], ese.robust=ese.robust, cp.robust=ci.properties[[2]][2],
                            width.robust=ci.properties[[1]][2], event.mean = event.mean, num.removed = num.removed, row.names = paste0("Strata=", strata))
  return(final.table)
}




## Run simulations ##
# Select cluster scenario
if (inputClusterType == 1){
  results.copula <- foreach(i=1:simCount, .packages=c('randomizr', 'copula')) %dopar% generateData.copula(n, K=3, inputTau, beta, randomization, strataStatus, censor.rate)
  clusterString <- "Small"
}else if (inputClusterType == 2){
  results.copula <- foreach(i=1:simCount, .packages=c('randomizr', 'copula')) %dopar% generateData.copula(n, K=300, inputTau, beta, randomization, strataStatus, censor.rate)
  clusterString <- "Large"
}else if (inputClusterType == 3){
  cnk <- c(rep(3, 5000), rep(300, 50)) # Custom nk vector for mixed clusters
  results.copula <- foreach(i=1:simCount, .packages=c('randomizr', 'copula')) %dopar% generateData.copula(n, NULL, inputTau, beta, randomization, strataStatus, censor.rate, custom_nk = cnk)
  clusterString <- "Mixed"
}

unstratifiedStats <- calcStats(results.copula, strata=FALSE, genCI)
stratifiedStats <- calcStats(results.copula, strata=TRUE, genCI)
final.table <- rbind(unstratifiedStats, stratifiedStats)


# Print out final table
write.csv(final.table, paste0("Tau", inputTau, "Cluster", clusterString, "Random", randomization, "Strata", strataStatus, ".csv"))




