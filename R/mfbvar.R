source("mat_gen_BMF.R")
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
  y_full <- load("data.dat")
  
  y_mod_full <- y_full[vec1,]
  data <-y_full[,1:(nobs + 1)]
  data_t <- t(data)
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

