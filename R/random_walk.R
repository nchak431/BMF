library(matrixcalc)
library(clusterGeneration)
library(fBasics)
library(LPTime)
library(StMoMo)
#==============================input======================================
k1 <- 3
k2 <- 30
nobs <- 100
horizon <- 20
#============================
# Load the data
y_full <- load("data.dat")

data <- y_full[,1:(nobs + 1)]
forecast_true_m1 <- y_full[,(nobs + 1 + 1):(nobs + 1 + horizon)]

for(i in 1:((3*k1)+k2)){
  data[i,] <- ts(data[i,],start = c(1990,1),frequency = 12)
}
#==============================
rw <- mrwd(data)
drift <- rw$drift

forecast <- list()
forecast[[1]] <- as.vector(drift) + data[,nobs+1] 
for ( t in 2:horizon)
{
  forecast[[t]] <- as.vector(drift) + forecast[[t-1]]
  
}
forecast <- do.call(cbind,forecast)
forecast <- as.matrix(forecast) 
err <- numeric()
for ( j in 1:horizon)
{
  err[j] <- sum((forecast_true_m1[,j][((3*k1)+1):((3*k1)+k2)]-forecast[,j][((3*k1)+1):((3*k1)+k2)])^2)/k2 #forecast error
}