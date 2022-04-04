source("mat_gen_BMF.R")

library(matrixcalc)
library(clusterGeneration)
library(fBasics)
library(midasr)

#==============================input======================================
k1 <- 3
k2 <- 30
nobs <- 101
horizon <- 20
#============================================================================
vec1 <- permutation(k1,k2)[[2]]

# Load the data
y_full <- load("data.dat")

data_mod_full <- y_full[vec1,]
data <-y_full[,1:(nobs + 1)]
data_t <- t(data)
data_mod <-as.matrix(y_mod_full)[,1:(nobs + 1)]
forecast_true_m1 <- y_full[,(nobs + 1 + 1):(nobs + 1 + horizon)]

mf_list <- list() 
for(i in 1:(k1+k2))  
{
  if(i <= k1) {
    ind <- c((3*i)-2,(3*i)-1,(3*i))  
    sub <- data_mod_full[ind,]
    mf_list[[i]] <- vec(sub)  
  }else {
    sub <- matrix(data_mod_full[((2*k1)+i),],ncol=1)
    mf_list[[i]] <- vec(sub)      
  } 
}
alfred_to_ts <- function(x,freq){
  ts(x[,1],frequency=freq)
}
mf <- mapply(alfred_to_ts,x=mf_list,freq=c(rep(12,k1),rep(4,k2)))

######################################################
# Data preparation

x1 <- window(mf[[1]], start = c(1, 1), end = c(26, 6))
x2 <- window(mf[[2]], start = c(1, 1), end = c(26, 6))
x3 <- window(mf[[3]], start = c(1, 1), end = c(26, 6))
x4 <- window(mf[[4]], start = c(1, 1), end = c(26, 2))
x5 <- window(mf[[5]], start = c(1, 1), end = c(26, 2))
x6 <- window(mf[[6]], start = c(1, 1), end = c(26, 2))
x7 <- window(mf[[7]], start = c(1, 1), end = c(26, 2))
x8 <- window(mf[[8]], start = c(1, 1), end = c(26, 2))
x9 <- window(mf[[9]], start = c(1, 1), end = c(26, 2))
x10 <- window(mf[[10]], start = c(1, 1), end = c(26, 2))
x11 <- window(mf[[11]], start = c(1, 1), end = c(26, 2))
x12 <- window(mf[[12]], start = c(1, 1), end = c(26, 2))
x13 <- window(mf[[13]], start = c(1, 1), end = c(26, 2))
x14 <- window(mf[[14]], start = c(1, 1), end = c(26, 2))
x15 <- window(mf[[15]], start = c(1, 1), end = c(26, 2))
x16 <- window(mf[[16]], start = c(1, 1), end = c(26, 2))
x17 <- window(mf[[17]], start = c(1, 1), end = c(26, 2))
x18 <- window(mf[[18]], start = c(1, 1), end = c(26, 2))
x19 <- window(mf[[19]], start = c(1, 1), end = c(26, 2))
x20 <- window(mf[[20]], start = c(1, 1), end = c(26, 2))
x21 <- window(mf[[21]], start = c(1, 1), end = c(26, 2))
x22 <- window(mf[[22]], start = c(1, 1), end = c(26, 2))
x23 <- window(mf[[23]], start = c(1, 1), end = c(26, 2))
x24 <- window(mf[[24]], start = c(1, 1), end = c(26, 2))
x25 <- window(mf[[25]], start = c(1, 1), end = c(26, 2))
x26 <- window(mf[[26]], start = c(1, 1), end = c(26, 2))
x27 <- window(mf[[27]], start = c(1, 1), end = c(26, 2))
x28 <- window(mf[[28]], start = c(1, 1), end = c(26, 2))
x29 <- window(mf[[29]], start = c(1, 1), end = c(26, 2))
x30 <- window(mf[[30]], start = c(1, 1), end = c(26, 2))
x31 <- window(mf[[31]], start = c(1, 1), end = c(26, 2))
x32 <- window(mf[[32]], start = c(1, 1), end = c(26, 2))
x33 <- window(mf[[33]], start = c(1, 1), end = c(26, 2))

fulldata <- list()
for(i in 1: k2){
  fulldata[[i]] <- list(x1 = window(mf[[1]], start = c(1, 1), end = c(28, 12)),x2 = window(mf[[2]], start = c(1, 1), end = c(28, 12)),
                        x3 = window(mf[[3]], start = c(1, 1), end = c(28, 12)),x4 = window(mf[[4]], start = c(1, 1), end = c(28, 4)),
                        x5 = window(mf[[5]], start = c(1, 1), end = c(28, 4)),x6= window(mf[[6]], start = c(1, 1), end = c(28, 4)),
                        x7 = window(mf[[7]], start = c(1, 1), end = c(28, 4)),x8 = window(mf[[8]], start = c(1, 1), end = c(28, 4)),
                        x9 = window(mf[[9]], start = c(1, 1), end = c(28, 4)),x10 = window(mf[[10]], start = c(1, 1), end = c(28, 4)),
                        x11 = window(mf[[11]], start = c(1, 1), end = c(28, 4)),x12 = window(mf[[12]], start = c(1, 1), end = c(28, 4)),
                        x13 = window(mf[[13]], start = c(1, 1), end = c(28, 4)),x14 = window(mf[[14]], start = c(1, 1), end = c(28, 4)),
                        x15 = window(mf[[15]], start = c(1, 1), end = c(28, 4)),x16 = window(mf[[16]], start = c(1, 1), end = c(28, 4)),
                        x17 = window(mf[[17]], start = c(1, 1), end = c(28, 4)),x18 = window(mf[[18]], start = c(1, 1), end = c(28, 4)),
                        x19 = window(mf[[19]], start = c(1, 1), end = c(28, 4)),x20 = window(mf[[20]], start = c(1, 1), end = c(28, 4)),
                        x21 = window(mf[[21]], start = c(1, 1), end = c(28, 4)),x22 = window(mf[[22]], start = c(1, 1), end = c(28, 4)),
                        x23 = window(mf[[23]], start = c(1, 1), end = c(28, 4)),x24 = window(mf[[24]], start = c(1, 1), end = c(28, 4)),
                        x25 = window(mf[[25]], start = c(1, 1), end = c(28, 4)),x26 = window(mf[[26]], start = c(1, 1), end = c(28, 4)),
                        x27 = window(mf[[27]], start = c(1, 1), end = c(28, 4)),x28 = window(mf[[28]], start = c(1, 1), end = c(28, 4)),
                        x29 = window(mf[[29]], start = c(1, 1), end = c(28, 4)),x30 = window(mf[[30]], start = c(1, 1), end = c(28, 4)),
                        x31 = window(mf[[31]], start = c(1, 1), end = c(28, 4)),x32 = window(mf[[32]], start = c(1, 1), end = c(28, 4)),
                        x33 = window(mf[[33]], start = c(1, 1), end = c(28, 4)),y = window(mf[[k1+i]], start = c(1, 1), end = c(28, 4)))
}

#==============================================================
# For forecasts corresponding to h=1,2,3,... 
midas_restr_s1 <- list()
midas_unrestr_s1 <- list()
forecast_restr_s1 <- matrix(0,nrow=k2,ncol=horizon)
forecast_unrestr_s1 <- matrix(0,nrow=k2,ncol=horizon)

for(i in 1:k2){
  y <- window(mf[[k1+i]], start = c(1, 1), end = c(26, 2))
  #Restricted MIDAS
  midas_restr_s1[[i]] <- midas_r(y ~ mls(x1,3:5,3,nealmon)+mls(x2,3:5,3,nealmon)+mls(x3,3:5,3,nealmon)+
                                   mls(x4,1,1,nealmon)+ mls(x5,1,1,nealmon)+mls(x6,1,1,nealmon)+
                                   mls(x7,1,1,nealmon)+ mls(x8,1,1,nealmon)+mls(x9,1,1,nealmon)+
                                   mls(x10,1,1,nealmon)+ mls(x11,1,1,nealmon)+mls(x12,1,1,nealmon)+
                                   mls(x13,1,1,nealmon)+ mls(x14,1,1,nealmon)+mls(x15,1,1,nealmon)+
                                   mls(x16,1,1,nealmon)+ mls(x17,1,1,nealmon)+mls(x18,1,1,nealmon)+
                                   mls(x19,1,1,nealmon)+ mls(x20,1,1,nealmon)+mls(x21,1,1,nealmon)+
                                   mls(x22,1,1,nealmon)+ mls(x23,1,1,nealmon)+mls(x24,1,1,nealmon)+
                                   mls(x25,1,1,nealmon)+ mls(x26,1,1,nealmon)+mls(x27,1,1,nealmon)+
                                   mls(x28,1,1,nealmon)+ mls(x29,1,1,nealmon)+mls(x30,1,1,nealmon)+
                                   mls(x31,1,1,nealmon)+ mls(x32,1,1,nealmon)+mls(x33,1,1,nealmon),
                                 start = list(x1= c(1,-0.5), x2 = c(1,-0.5),x3 = c(1,-0.5),x4 = c(1,-0.5),
                                              x5 = c(1,-0.5),x6 = c(1,-0.5),x7 = c(1,-0.5),x8 = c(1,-0.5),
                                              x9 = c(1,-0.5),x10 = c(1,-0.5),x11 = c(1,-0.5),x12 = c(1,-0.5),
                                              x13 = c(1,-0.5),x14 = c(1,-0.5),x15 = c(1,-0.5),x16 = c(1,-0.5),
                                              x17 = c(1,-0.5),x18 = c(1,-0.5),x19 = c(1,-0.5),x20 = c(1,-0.5),
                                              x21 = c(1,-0.5),x22 = c(1,-0.5),x23 = c(1,-0.5),x24 = c(1,-0.5),
                                              x25 = c(1,-0.5),x26 = c(1,-0.5),x27 = c(1,-0.5),x28 = c(1,-0.5),
                                              x29 = c(1,-0.5),x30 = c(1,-0.5),x31 = c(1,-0.5),x32 = c(1,-0.5),x33 = c(1,-0.5)), Ofunction = "optim", method = "Nelder-Mead")
  #Unrestricted MIDAS
  midas_unrestr_s1[[i]] <- midas_r(y ~ mls(x1,3:5,3)+mls(x2,3:5,3)+mls(x3,3:5,3)+
                                     mls(x4,1,1)+ mls(x5,1,1)+mls(x6,1,1)+
                                     mls(x7,1,1)+ mls(x8,1,1)+mls(x9,1,1)+
                                     mls(x10,1,1)+ mls(x11,1,1)+mls(x12,1,1)+
                                     mls(x13,1,1)+ mls(x14,1,1)+mls(x15,1,1)+
                                     mls(x16,1,1)+ mls(x17,1,1)+mls(x18,1,1)+
                                     mls(x19,1,1)+ mls(x20,1,1)+mls(x21,1,1)+
                                     mls(x22,1,1)+ mls(x23,1,1)+mls(x24,1,1)+
                                     mls(x25,1,1)+ mls(x26,1,1)+mls(x27,1,1)+
                                     mls(x28,1,1)+ mls(x29,1,1)+mls(x30,1,1)+
                                     mls(x31,1,1)+ mls(x32,1,1)+mls(x33,1,1),
                                   start = NULL)
  
  
  full <- fulldata[[i]]
  insample <- 1:length(y)
  outsample <- (1:length(full$y))[-insample]
  
  avgf_restr_s1 <- average_forecast(list(midas_restr_s1[[i]]), data = full, insample = insample, outsample = outsample)
  avgf_unrestr_s1 <- average_forecast(list(midas_unrestr_s1[[i]]), data = full, insample = insample, outsample = outsample)
  
  forecast_restr_s1[i,] <- avgf_restr_s1$forecast
  forecast_unrestr_s1[i,] <- avgf_unrestr_s1$forecast
}

err_low_restr_s1 <- numeric()
err_low_unrestr_s1 <- numeric()
for ( j in 1:horizon)
{
  err_low_restr_s1[j] <- sum((forecast_true[,j]-forecast_restr_s1[,j])^2)/k2 #Forecast error for restricted MIDAS
  err_low_unrestr_s1[j] <- sum((forecast_true[,j]-forecast_unrestr_s1[,j])^2)/k2 #Forecast error for unrestricted MIDAS
}
##############################################################################
# For forecasts corresponding to h=2/3,5/3,... 
midas_restr_m1 <- list()
midas_unrestr_m1 <- list()
forecast_restr_m1 <- matrix(0,nrow=k2,ncol=horizon)
forecast_unrestr_m1 <- matrix(0,nrow=k2,ncol=horizon)

for(i in 1:k2){
  y <- window(mf[[k1+i]], start = c(1, 1), end = c(26, 2))
  midas_restr_m1[[i]] <- midas_r(y ~ mls(x1,2:4,3,nealmon)+mls(x2,2:4,3,nealmon)+mls(x3,2:4,3,nealmon)+
                                   mls(x4,1,1,nealmon)+ mls(x5,1,1,nealmon)+mls(x6,1,1,nealmon)+
                                   mls(x7,1,1,nealmon)+ mls(x8,1,1,nealmon)+mls(x9,1,1,nealmon)+
                                   mls(x10,1,1,nealmon)+ mls(x11,1,1,nealmon)+mls(x12,1,1,nealmon)+
                                   mls(x13,1,1,nealmon)+ mls(x14,1,1,nealmon)+mls(x15,1,1,nealmon)+
                                   mls(x16,1,1,nealmon)+ mls(x17,1,1,nealmon)+mls(x18,1,1,nealmon)+
                                   mls(x19,1,1,nealmon)+ mls(x20,1,1,nealmon)+mls(x21,1,1,nealmon)+
                                   mls(x22,1,1,nealmon)+ mls(x23,1,1,nealmon)+mls(x24,1,1,nealmon)+
                                   mls(x25,1,1,nealmon)+ mls(x26,1,1,nealmon)+mls(x27,1,1,nealmon)+
                                   mls(x28,1,1,nealmon)+ mls(x29,1,1,nealmon)+mls(x30,1,1,nealmon)+
                                   mls(x31,1,1,nealmon)+ mls(x32,1,1,nealmon)+mls(x33,1,1,nealmon),
                                 start = list(x1 = c(1,-0.5),x2 = c(1,-0.5),x3 = c(1,-0.5),x4 = c(1,-0.5),
                                              x5 = c(1,-0.5),x6 = c(1,-0.5),x7 = c(1,-0.5),x8 = c(1,-0.5),
                                              x9 = c(1,-0.5),x10 = c(1,-0.5),x11 = c(1,-0.5),x12 = c(1,-0.5),
                                              x13 = c(1,-0.5),x14 = c(1,-0.5),x15 = c(1,-0.5),x16 = c(1,-0.5),
                                              x17 = c(1,-0.5),x18 = c(1,-0.5),x19 = c(1,-0.5),x20 = c(1,-0.5),
                                              x21 = c(1,-0.5),x22 = c(1,-0.5),x23 = c(1,-0.5),x24 = c(1,-0.5),
                                              x25 = c(1,-0.5),x26 = c(1,-0.5),x27 = c(1,-0.5),x28 = c(1,-0.5),
                                              x29 = c(1,-0.5),x30 = c(1,-0.5),x31 = c(1,-0.5),x32 = c(1,-0.5),
                                              x33 = c(1,-0.5)), Ofunction = "optim", method = "Nelder-Mead")
  midas_unrestr_m1[[i]] <- midas_r(y ~ mls(x1,2:4,3)+mls(x2,2:4,3)+mls(x3,2:4,3)+
                                     mls(x4,1,1)+ mls(x5,1,1)+mls(x6,1,1)+
                                     mls(x7,1,1)+ mls(x8,1,1)+mls(x9,1,1)+
                                     mls(x10,1,1)+ mls(x11,1,1)+mls(x12,1,1)+
                                     mls(x13,1,1)+ mls(x14,1,1)+mls(x15,1,1)+
                                     mls(x16,1,1)+ mls(x17,1,1)+mls(x18,1,1)+
                                     mls(x19,1,1)+ mls(x20,1,1)+mls(x21,1,1)+
                                     mls(x22,1,1)+ mls(x23,1,1)+mls(x24,1,1)+
                                     mls(x25,1,1)+ mls(x26,1,1)+mls(x27,1,1)+
                                     mls(x28,1,1)+ mls(x29,1,1)+mls(x30,1,1)+
                                     mls(x31,1,1)+ mls(x32,1,1)+mls(x33,1,1),
                                   start = NULL)
  
  
  full <- fulldata[[i]]
  insample <- 1:length(y)
  outsample <- (1:length(full$y))[-insample]
  
  avgf_restr_m1 <- average_forecast(list(midas_restr_m1[[i]]), data = full, insample = insample, outsample = outsample)
  avgf_unrestr_m1 <- average_forecast(list(midas_unrestr_m1[[i]]), data = full, insample = insample, outsample = outsample)
  
  forecast_restr_m1[i,] <- avgf_restr_m1$forecast
  forecast_unrestr_m1[i,] <- avgf_unrestr_m1$forecast
}

err_low_restr_m1 <- numeric()
err_low_unrestr_m1 <- numeric()
for ( j in 1:horizon)
{
  err_low_restr_m1[j] <- sum((forecast_true[,j]-forecast_restr_m1[,j])^2)/k2
  err_low_unrestr_m1[j] <- sum((forecast_true[,j]-forecast_unrestr_m1[,j])^2)/k2
}


##############################################################################
# For forecasts corresponding to h=1/3,4/3,... 

midas_restr_m2 <- list()
midas_unrestr_m2 <- list()
forecast_restr_m2 <- matrix(0,nrow=k2,ncol=horizon)
forecast_unrestr_m2 <- matrix(0,nrow=k2,ncol=horizon)

for(i in 1:k2){
  y <- window(mf[[k1+i]], start = c(1, 1), end = c(26, 2))
  midas_restr_m2[[i]] <- midas_r(y ~ mls(x1,1:3,3,nealmon)+mls(x2,1:3,3,nealmon)+mls(x3,1:3,3,nealmon)+
                                   mls(x4,1,1,nealmon)+ mls(x5,1,1,nealmon)+mls(x6,1,1,nealmon)+
                                   mls(x7,1,1,nealmon)+ mls(x8,1,1,nealmon)+mls(x9,1,1,nealmon)+
                                   mls(x10,1,1,nealmon)+ mls(x11,1,1,nealmon)+mls(x12,1,1,nealmon)+
                                   mls(x13,1,1,nealmon)+ mls(x14,1,1,nealmon)+mls(x15,1,1,nealmon)+
                                   mls(x16,1,1,nealmon)+ mls(x17,1,1,nealmon)+mls(x18,1,1,nealmon)+
                                   mls(x19,1,1,nealmon)+ mls(x20,1,1,nealmon)+mls(x21,1,1,nealmon)+
                                   mls(x22,1,1,nealmon)+ mls(x23,1,1,nealmon)+mls(x24,1,1,nealmon)+
                                   mls(x25,1,1,nealmon)+ mls(x26,1,1,nealmon)+mls(x27,1,1,nealmon)+
                                   mls(x28,1,1,nealmon)+ mls(x29,1,1,nealmon)+mls(x30,1,1,nealmon)+
                                   mls(x31,1,1,nealmon)+ mls(x32,1,1,nealmon)+mls(x33,1,1,nealmon),
                                 start = list(x1 = c(1,-0.5),x2 = c(1,-0.5),x3 = c(1,-0.5),x4 = c(1,-0.5),
                                              x5 = c(1,-0.5),x6 = c(1,-0.5),x7 = c(1,-0.5),x8 = c(1,-0.5),
                                              x9 = c(1,-0.5),x10 = c(1,-0.5),x11 = c(1,-0.5),x12 = c(1,-0.5),
                                              x13 = c(1,-0.5),x14 = c(1,-0.5),x15 = c(1,-0.5),x16 = c(1,-0.5),
                                              x17 = c(1,-0.5),x18 = c(1,-0.5),x19 = c(1,-0.5),x20 = c(1,-0.5),
                                              x21 = c(1,-0.5),x22 = c(1,-0.5),x23 = c(1,-0.5),x24 = c(1,-0.5),
                                              x25 = c(1,-0.5),x26 = c(1,-0.5),x27 = c(1,-0.5),x28 = c(1,-0.5),
                                              x29 = c(1,-0.5),x30 = c(1,-0.5),x31 = c(1,-0.5),x32 = c(1,-0.5),
                                              x33 = c(1,-0.5)), Ofunction = "optim", method = "Nelder-Mead")
  midas_unrestr_m2[[i]] <- midas_r(y ~ mls(x1,1:3,3)+mls(x2,1:3,3)+mls(x3,1:3,3)+
                                     mls(x4,1,1)+ mls(x5,1,1)+mls(x6,1,1)+
                                     mls(x7,1,1)+ mls(x8,1,1)+mls(x9,1,1)+
                                     mls(x10,1,1)+ mls(x11,1,1)+mls(x12,1,1)+
                                     mls(x13,1,1)+ mls(x14,1,1)+mls(x15,1,1)+
                                     mls(x16,1,1)+ mls(x17,1,1)+mls(x18,1,1)+
                                     mls(x19,1,1)+ mls(x20,1,1)+mls(x21,1,1)+
                                     mls(x22,1,1)+ mls(x23,1,1)+mls(x24,1,1)+
                                     mls(x25,1,1)+ mls(x26,1,1)+mls(x27,1,1)+
                                     mls(x28,1,1)+ mls(x29,1,1)+mls(x30,1,1)+
                                     mls(x31,1,1)+ mls(x32,1,1)+mls(x33,1,1),
                                   start = NULL)
  
  
  full <- fulldata[[i]]
  insample <- 1:length(y)
  outsample <- (1:length(full$y))[-insample]
  
  avgf_restr_m2 <- average_forecast(list(midas_restr_m2[[i]]), data = full, insample = insample, outsample = outsample)
  avgf_unrestr_m2 <- average_forecast(list(midas_unrestr_m2[[i]]), data = full, insample = insample, outsample = outsample)
  
  forecast_restr_m2[i,] <- avgf_restr_m2$forecast
  forecast_unrestr_m2[i,] <- avgf_unrestr_m2$forecast
}

err_low_restr_m2 <- numeric()
err_low_unrestr_m2 <- numeric()
for ( j in 1:horizon)
{
  err_low_restr_m2[j] <- sum((forecast_true[,j]-forecast_restr_m2[,j])^2)/k2
  err_low_unrestr_m2[j] <- sum((forecast_true[,j]-forecast_unrestr_m2[,j])^2)/k2
}
