source("data_gen.R")
source("functions.R")

library(matrixcalc)
library(clusterGeneration)
library(fBasics)
library(glmnet)
library(psych)
library(MASS)
library(Matrix)
library(glmnet)
library(mfbvar)
library(scoringRules)
#==============================input======================================
k1 <- 3
k2 <- 30
nobs <- 100
horizon <- 10
repl <- 1 #For a single replicate
#============================================================================
vec1 <- permutation(k1,k2)[[2]]

p <- (3*k1) + k2
p1 <- (k1+k2)
pred1 <- array(0, c(((3*(k1+k2))*horizon),3, repl))
pred2 <- array(0, c(((3*(k1+k2))*horizon),3, repl))
pred3 <- array(0, c(((3*(k1+k2))*horizon),3, repl))

crps_minn <- array(0, c(k2,(3*horizon), repl))
crps_ss <- array(0, c(k2,(3*horizon), repl))
crps_ssng <- array(0, c(k2,(3*horizon), repl))

logs_minn <- array(0, c(k2,(3*horizon), repl))
logs_ss <- array(0, c(k2,(3*horizon), repl))
logs_ssng <- array(0, c(k2,(3*horizon), repl))

ferr_mfb_low1 <- array(0, c((3*horizon),1, repl))
ferr_mfb_low2 <- array(0, c((3*horizon),1, repl))
ferr_mfb_low3 <- array(0, c((3*horizon),1, repl))

#For a single replicate
r=1
# Load the data
load(paste("y_full_test-", index, ".dat", sep=''))
y_full <- y_full
y_mod_full <- y_full[vec1,]
data <-y_full[,1:(nobs + 1)]
data_mod <-as.matrix(y_mod_full)[,1:(nobs + 1)]
forecast_true_m1 <- y_full[,(nobs + 1 + 1):(nobs + 1 + horizon)]
  
  arr <- mfb_rearr(k1,k2,horizon)
  mfb_compare <- mfb(data_mod,k1,k2,horizon,arr,forecast_true_m1)
  #Under Minnesota prior 
  pred1[,,r] <- as.matrix(mfb_compare$pred_1)
  crps_minn[,,r] <- mfb_compare$crps_minn
  logs_minn[,,r] <- mfb_compare$logs_minn
  ferr_mfb_low1[,,r] <- mfb_compare$ferr_mfb_low1
  #Under Steady-state prior
  pred2[,,r] <- as.matrix(mfb_compare$pred_2)
  crps_ss[,,r] <- mfb_compare$crps_ss
  logs_ss[,,r] <- mfb_compare$logs_ss
  ferr_mfb_low2[,,r] <- mfb_compare$ferr_mfb_low2
  #Under Hierarchical steady-state prior
  pred3[,,r] <- as.matrix(mfb_compare$pred_3)
  crps_ssng[,,r] <- mfb_compare$crps_ssng
  logs_ssng[,,r] <- mfb_compare$logs_ssng
  ferr_mfb_low3[,,r] <- mfb_compare$ferr_mfb_low3 
index <- 1 # Fix an index no. for each replicate.
save(crps_minn, file=paste("crps_minn_test-", index, ".dat", sep=''))
save(logs_minn, file=paste("logs_minn_test-", index, ".dat", sep=''))
save(crps_ss, file=paste("crps_ss_test-", index, ".dat", sep=''))
save(logs_ss, file=paste("logs_ss_test-", index, ".dat", sep=''))
save(crps_ssng, file=paste("crps_ssng_test-", index, ".dat", sep=''))
save(logs_ssng, file=paste("logs_ssng_test-", index, ".dat", sep=''))

save(ferr_mfb_low1, file=paste("ferr_low_minn_test-", index, ".dat", sep=''))
save(ferr_mfb_low2, file=paste("ferr_low_ss_test-", index, ".dat", sep=''))
save(ferr_mfb_low3, file=paste("ferr_low_ssng_test-", index, ".dat", sep=''))

# Repeat this for 50(or more) replicates and then combine the results.

# Final RMSE
repl <- 50
val_minn <- NULL
val_ss <- NULL
val_ssng <- NULL

for(i in 1:repl){
    load(paste("ferr_low_minn_test-", i, ".dat", sep=''))
    load(paste("ferr_low_ss_test-", i, ".dat", sep=''))
    load(paste("ferr_low_ssng_test-", i, ".dat", sep=''))
    err_minn <-ferr_mfb_low1[,,1]
    err_ss <-ferr_mfb_low2[,,1]
    err_ssng <-ferr_mfb_low3[,,1]
    val_minn <- rbind(val_minn,err_minn) 
    val_ss <- rbind(val_ss,err_ss) 
    val_ssng <- rbind(val_ssng,err_ssng) 
}
minn <- apply(val_minn,2,mean)
ss <- apply(val_ss,2,mean)
ssng <- apply(val_ssng,2,mean)

RMSE_minn <- sqrt(minn)
RMSE_ss <- sqrt(ss)
RMSE_ssng <- sqrt(ssng)

#Final CRPS and LPS

avg_crps_minn <- NULL
avg_logs_minn <- NULL
avg_crps_ss <- NULL
avg_logs_ss <- NULL
avg_crps_ssng <- NULL
avg_logs_ssng <- NULL
for(ind in 1:repl){
  load(paste("crps_minn_test-", ind, ".dat", sep=''))               
  load(paste("logs_minn_test-", ind, ".dat", sep=''))
  load(paste("crps_ss_test-", ind, ".dat", sep=''))               
  load(paste("logs_ss_test-", ind, ".dat", sep=''))
  
  crps_minn <- crps_minn[,,1]
  logs_minn <- logs_minn[,,1]
  crps_ss <- crps_ss[,,1]
  logs_ss <- logs_ss[,,1]
  
  minn_crps <- numeric()
  minn_logs <- numeric()
  for ( j in 1:(3*horizon))
  {
    minn_crps[j] <- mean(crps_minn[,j])
    minn_logs[j] <- mean(logs_minn[,j])
  } 
  avg_crps_minn <- rbind(avg_crps_minn,minn_crps)
  avg_logs_minn <- rbind(avg_logs_minn,minn_logs)
  
  ss_crps <- numeric()
  ss_logs <- numeric()
  for ( j in 1:(3*horizon))
  {
    ss_crps[j] <- mean(crps_ss[,j])
    ss_logs[j] <- mean(logs_ss[,j])
  } 
  avg_crps_ss <- rbind(avg_crps_ss,ss_crps)
  avg_logs_ss <- rbind(avg_logs_ss,ss_logs)
  
 
}

final_crps_minn <- apply(avg_crps_minn,2,mean)
final_logs_minn <- apply(avg_logs_minn,2,mean)

final_crps_ss <- apply(avg_crps_ss,2,mean)
final_logs_ss <- apply(avg_logs_ss,2,mean)


for(ind in 1:repl){
   load(paste("crps_ssng_test-", ind, ".dat", sep=''))               
  load(paste("logs_ssng_test-", ind, ".dat", sep=''))
  crps_ssng <- crps_ssng[,,1]
  logs_ssng <- logs_ssng[,,1]
  ssng_crps <- numeric()
  ssng_logs <- numeric()
  for ( j in 1:(3*horizon))
  {
    ssng_crps[j] <- mean(crps_ssng[,j])
    ssng_logs[j] <- mean(logs_ssng[,j])
  } 
  avg_crps_ssng <- rbind(avg_crps_ssng,ssng_crps)
  avg_logs_ssng <- rbind(avg_logs_ssng,ssng_logs)
  
}

final_crps_ssng <- apply(avg_crps_ssng,2,mean)
final_logs_ssng <- apply(avg_logs_ssng,2,mean)


crps <- cbind(final_crps_minn,final_crps_ss,final_crps_ssng)
logs <- cbind(final_logs_minn,final_logs_ss,final_logs_ssng)

