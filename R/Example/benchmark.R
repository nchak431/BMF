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
load(paste("y_full_test-", index, ".dat", sep=''))
y_full <- y_full
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
index <- 1 # Fix an index no. for each replicate.
save(err, file=paste("err_test-", index, ".dat", sep=''))

# Repeat this for 50(or more) replicates and then combine the results.

# Final RMSE
repl <- 50

m <- array(0, c((horizon),1, repl))
for(i in 1:repl){
    load(paste("err_th0.2_test-", i, ".dat", sep=''))
    
    m[,,i]<-  err     
}
final <- apply(m,c(1,2),mean)
RMSE <- sqrt(final)
