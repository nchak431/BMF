#Prepare Data 2

library(dplyr)
library(alfred)
library(matrixcalc)
library(MASS)
library(Matrix)
library(glmnet)

variables <- c("GDPC1","PCECC96","GPDIC1","GCEC1","IPDBS","RCPHBS",
               "OPHPBS","TABSHNO","TLBSNNCB",
               "CPIAUCSL","CPIULFSL","CPIAPPSL","PCEPI","CPITRNSL","DSERRG3M086SBEA","MZMSL","TOTRESNS","BUSLOANS","NONREVSL",
               "PPICMM","DTCTHFNM","CES3000000008","CES0600000008","CUSR0000SA0L2",
               "DTCOLNVHFNM","DNDGRG3M086SBEA","CPIMEDSL",#tcode 6
               "INDPRO","IPDMAT","RPI","IPBUSEQ","W875RX1","DPCERA3M086SBEA","USTPU",	"CE16OV",	"CLF16OV",
               "IPNMAT","IPFINAL","IPFPNSS","IPFUELS","IPNCONGD","USFIRE",
               "USCONS","USGOOD","M2REAL","DMANEMP",#tcode 5
               "FEDFUNDS","TB3MS","gs1","gs5","gs10","TB6MS","UNRATE","AAA","BAA","UEMPMEAN","AWOTMAN") #tcode 2
              
out <- list()
out <- lapply(variables,get_alfred_series,observation_start="1959-04-01",
              observation_end="2019-09-01",realtime_start="2021-09-10",realtime_end="2021-09-10")#Q2 1959 to Q3 2019

#"SP500"
sp <-read.csv("S&P_v1.csv")
sp_dataframe <- data.frame(date=rnorm(length(sp)),realtime_period=rnorm(length(sp)),sp=sp)

out[[l+1]] <- sp_dataframe

#Transformations

new <- list()
set1 <- c(1:4,6:9)
for(i in set1){
  temp <- as.numeric(out[[i]][,3])
  diff <- 100*diff(log(temp))
  new[[i]] <- as.data.frame(diff[-c(1:2)])
}

set2 <- c(5)
for(i in set2){
  temp <- as.numeric((out[[i]][,3]))
  diff <-  100* diff(log(temp),diff=2)
  new[[i]] <- as.data.frame(diff[-1])
}
variables[47]

set3 <- c(10:27)
for(i in set3){
  temp <- as.numeric((out[[i]][,3]))
  diff <-  100* diff(log(temp),diff=2)
  new[[i]] <- as.data.frame(diff[-c(1:7)])
}
set4 <- c(28:46)
for(i in set4){
  temp <- as.numeric(out[[i]][,3])
  diff <- 100*diff(log(temp))
  new[[i]] <- as.data.frame(diff[-c(1:8)])
}

set5 <- c(47:57)
for(i in set5){
  temp <- as.numeric((out[[i]][,3]))
  diff <- diff(temp)
  new[[i]] <- as.data.frame(diff[-c(1:8)])
}


set6 <- c((l+1))
for(i in set6){
  temp <- as.numeric((out[[i]][,3]))
  diff <- 100*diff(log(temp))
  new[[i]] <- as.data.frame(diff[-c(1:8)])
}

k1 <- 49
k2 <- 9
new1 <- list()
for(i in 1:k1){
  new1[[i]] <- new[[k2+i]]  
}
for(i in (k1+1):(k1+k2)){
  new1[[i]] <- new[[i-k1]]  
}

variables1 <- c("CPIAUCSL","CPIULFSL","CPIAPPSL","PCEPI","CPITRNSL","DSERRG3M086SBEA","MZMSL","TOTRESNS","BUSLOANS","NONREVSL",
                "PPICMM","DTCTHFNM","CES3000000008","CES0600000008","CUSR0000SA0L2",
                "DTCOLNVHFNM","DNDGRG3M086SBEA","CPIMEDSL",#tcode 6
                "INDPRO","IPDMAT","RPI","IPBUSEQ","W875RX1","DPCERA3M086SBEA","USTPU",	"CE16OV",	"CLF16OV",
                "IPNMAT","IPFINAL","IPFPNSS","IPFUELS","IPNCONGD","USFIRE",
                "USCONS","USGOOD","M2REAL","DMANEMP",#tcode 5
                "FEDFUNDS","TB3MS","gs1","gs5","gs10","TB6MS","UNRATE","AAA","BAA","UEMPMEAN","AWOTMAN", #tcode 2
               "sp", #tcode 5
               "GDPC1","PCECC96","GPDIC1","GCEC1","IPDBS","RCPHBS",
               "OPHPBS","TABSHNO","TLBSNNCB") 


length(variables1)
names(new1) <- variables1

alfred_to_ts <- function(x,freq){
  ts(x[,1],start=c(1960,1),frequency=freq)
}
mf_list <- list()
for(i in 1:k1){
  mf_list[[i]]<- ts(new1[[i]][,1],start=c(1960,1),frequency=12) 
}
for(i in (k1+1):(k1+k2) ){
  mf_list[[i]]<- ts(new1[[i]][,1],start=c(1960,1),frequency=4) 
}


s1 <- NULL
for(i in 1:length(mf_list[[(k1+1)]])){
  m <- c((3*(i-1))+1)  
  s1 <-c(s1,m)
}

s2 <- NULL
for(i in 1:length(mf_list[[(k1+1)]])){
  m <- c((3*(i-1))+2)  
  s2 <-c(s2,m)
}
s3 <- NULL
for(i in 1:length(mf_list[[(k1+1)]])){
  m <- c((3*(i-1))+3)  
  s3 <-c(s3,m)
}
hf1 <- NULL
hf2 <- NULL
hf3 <- NULL
for(i in 1:k1){
  hf3 <- rbind(hf3, mf_list[[i]][s3])
  hf1 <- rbind(hf1, mf_list[[i]][s1])
  hf2 <- rbind(hf2, mf_list[[i]][s2])
}
hf <- rbind(hf3,hf1,hf2)
lf <- NULL
for(i in (k1+1):(k1+k2)){
  lf <- rbind(lf, mf_list[[i]])
}
data_full <- rbind(hf,lf)
write.csv(data_full,"data2.csv")