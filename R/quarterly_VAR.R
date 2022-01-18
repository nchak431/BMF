library(LPTime)
data <- read.csv("data1.csv") # Read real data (Data 1), would be similar for Data 2/Data 3
data <- data[,-1]
k1 <- 24
k2<- 9
data_hh <- data[1:(3*k1),]
data_ll <- data[((3*k1)+1):((3*k1)+k2),]

# Prepare quarterly data with equal weights
m1 <- NULL
for(i in 1:k1){
  m = c(i,(k1+i),((2*k1)+i))
  m1 <- c(m1,m)
}
data_hh <- data_hh[m1,] 

mat <- matrix(0,nrow=k1,ncol=238)
for(i in 1:238){
  m <- NULL
  for(j in 1:k1){
    temp <- sum(data_hh[,i][((3*j)-2)],data_hh[,i][((3*j)-1)],data_hh[,i][(3*j)])/3
    m <- c(m,temp)
  }
  mat[,i] <- m
  
}
data_full <- rbind(mat,data_ll)
#=====================================================

index=0
horizon <- 10
num <- 50 # no. of estimation samples
ntrain <- 228 # no. of samples in train data set (Jan 1960 - Dec 2016)

for(k in 1:num){
  index <- index + 1
  nobs <- ntrain - index
  for(i in 1: (k1+k2))
  {
    data_full[i,] <- data_full[i,] - mean(data_full[i,]) 
    
  }
  data <-as.matrix(data_full[,1:(nobs + 1)])
  forecast_true_m1 <- data_full[,(nobs + 1 + 1):(nobs + 1 + horizon)]
  data <- t(data)
  for(i in 1:(k1+k2)){
    data[,i] <- ts(data[,i],start = c(1990,1),frequency = 12)
  }
  data
  est <- VAR(data,p=1)
  est_mat <- as.matrix(est)
  
  forecast <- list()
  forecast[[1]] <- as.vector((est_mat %*% data_full[,nobs+1]))# + err)
  for ( t in 2:horizon)
  {
    forecast[[t]] <- as.vector((est_mat %*% forecast[[t-1]]))# + err)
    
  }
  forecast <- do.call(cbind,forecast)
  forecast <- as.matrix(forecast) 
  err_low <- numeric() # forecast error
  for ( j in 1:horizon)
  {
    err_low[j] <- sum((forecast_true_m1[,j][((k1)+1):((k1)+k2)]-forecast[,j][((k1)+1):((k1)+k2)])^2)/k2
  }
  save(err_low, file=paste("err_low-", index, ".dat", sep=''))
}
#COMBINING 5O ERRORS AND CALCULATE RMSE
m <- array(0, c(horizon,1, repl))
for(i in 1:num){
  load(paste("err_low-", i, ".dat", sep=''))
  m[,,i]<-  err_low     
}
final <- apply(m,c(1,2),mean)
RMSE <- sqrt(final)
